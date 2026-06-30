#!/usr/bin/env python3
#
#  pocket_obs2rnx.py - Convert pocket.log $OBS lines to RINEX 3.03 OBS file.
#
#  The resulting RINEX file can be used as a rover observation file for RTK
#  post-processing with rnx2rtkp (RTKLIB) and a base RINEX from a CORS station.
#
#  Usage:
#    python3 pocket_obs2rnx.py pocket.log [-o rover.obs] [options]
#
#  Options:
#    -o FILE       output RINEX 3 OBS file (default: same stem as log + .obs)
#    --sys G,R,E,J,C  satellite systems to include (default: all found)
#    --marker NAME    marker name for RINEX header (default: ROVER)
#
#  Example (drone experiment):
#    python3 pocket_obs2rnx.py work/20260701.../pocket.log -o rover.obs
#    # Download base RINEX from GA CORS (e.g. CANB):
#    #   https://gnss.ga.gov.au/rinex/  → CANB, select date + hour
#    /path/to/rnx2rtkp -k rtk_kinematic.conf rover.obs base.obs base.nav -o rtk.pos
#    python3 python/pocket_pos_plot.py rtk.pos --ref LAT LON HGT
#

import sys, argparse, os, math
from collections import defaultdict

# ── constants ──────────────────────────────────────────────────────────────────

RINEX_VER = '3.03'
PROG_NAME = 'pocket_obs2rnx.py'

_WGS84_A  = 6378137.0
_WGS84_E2 = 2 / 298.257223563 - (1 / 298.257223563) ** 2

# Map first char of sat ID → RINEX 3 system code
SYS_CHAR = {'G': 'G', 'R': 'R', 'E': 'E', 'J': 'J', 'C': 'C', 'S': 'S'}

# Preferred obs-type order within a system (only types present in data appear)
OBS_ORDER = ['C1C','L1C','D1C','S1C',
             'C1P','L1P','D1P','S1P',
             'C1S','L1S','D1S','S1S',
             'C1L','L1L','D1L','S1L',
             'C1E','L1E','D1E','S1E',
             'C2S','L2S','D2S','S2S',
             'C2L','L2L','D2L','S2L',
             'C2P','L2P','D2P','S2P',
             'C2C','L2C','D2C','S2C',
             'C5Q','L5Q','D5Q','S5Q',
             'C5P','L5P','D5P','S5P',
             'C5I','L5I','D5I','S5I',
             'C6S','L6S','D6S','S6S',
             'C6L','L6L','D6L','S6L']

# ── parsing ────────────────────────────────────────────────────────────────────

def parse_log(path):
    """Parse $OBS and first $POS from pocket.log.

    Returns:
        obs_records: list of dicts with keys t, epoch, sat, code, P, L, D, SNR, LLI
        pos_xyz:     approximate ECEF (from first PPP $POS), or None
    """
    obs_records = []
    pos_xyz = None

    with open(path, errors='replace') as f:
        for line in f:
            if line.startswith('$OBS,'):
                p = line.strip().split(',')
                if len(p) < 15:
                    continue
                try:
                    obs_records.append({
                        't':    float(p[1]),
                        'epoch': (int(float(p[2])), int(float(p[3])),
                                  int(float(p[4])), int(float(p[5])),
                                  int(float(p[6])), float(p[7])),
                        'sat':  p[8],
                        'code': p[9],
                        'SNR':  float(p[10]),
                        'P':    float(p[11]),
                        'L':    float(p[12]),
                        'D':    float(p[13]),
                        'LLI':  int(p[14]),
                    })
                except (ValueError, IndexError):
                    continue
            elif line.startswith('$POS,') and pos_xyz is None:
                p = line.strip().split(',')
                if len(p) >= 12 and int(p[11]) in (5, 6):
                    try:
                        lat, lon, hgt = float(p[8]), float(p[9]), float(p[10])
                        pos_xyz = llh_to_ecef(lat, lon, hgt)
                    except (ValueError, IndexError):
                        pass

    return obs_records, pos_xyz


def llh_to_ecef(lat_deg, lon_deg, hgt_m):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    N = _WGS84_A / math.sqrt(1.0 - _WGS84_E2 * math.sin(lat)**2)
    X = (N + hgt_m) * math.cos(lat) * math.cos(lon)
    Y = (N + hgt_m) * math.cos(lat) * math.sin(lon)
    Z = (N * (1.0 - _WGS84_E2) + hgt_m) * math.sin(lat)
    return X, Y, Z


# ── RINEX 3 formatting ─────────────────────────────────────────────────────────

def obs_code_to_rinex3(code2ch):
    """Map $OBS 2-char code → 4 RINEX-3 obs-type strings (C, L, D, S)."""
    return [f'C{code2ch}', f'L{code2ch}', f'D{code2ch}', f'S{code2ch}']


def snr_to_indicator(snr_dbhz):
    """Convert SNR in dB-Hz to RINEX signal-strength indicator (1-9, 0=unknown)."""
    if snr_dbhz <= 0:
        return 0
    # RINEX spec: 1=minimum, 9=maximum; 6 ≈ ≥ 42 dBHz
    v = int(snr_dbhz / 6)
    return min(max(v, 1), 9)


def fmt_obs(val, lli=0, snr_ind=0, is_phase=False):
    """Format a single 16-char RINEX 3 observable field."""
    if val == 0.0:
        return ' ' * 16
    s = f'{val:14.3f}'
    if is_phase:
        s += f'{lli & 3:1d}{snr_ind:1d}'
    else:
        s += '  '
    return s


# ── epoch grouping ─────────────────────────────────────────────────────────────

def group_epochs(obs_records, sys_filter):
    """Group obs_records by epoch and satellite.

    Returns:
        epochs:    list of (epoch_key, sat_dict) sorted by time
                   sat_dict: {sat → {rinex3_obs_type → (value, lli, snr_ind)}}
        all_types: {sys_char → sorted list of RINEX 3 obs types present}
    """
    # epoch_key → sat → code2ch → measurements
    raw = defaultdict(lambda: defaultdict(dict))

    for r in obs_records:
        sys_ch = r['sat'][0] if r['sat'] else ''
        if sys_filter and sys_ch not in sys_filter:
            continue
        if sys_ch not in SYS_CHAR:
            continue
        if r['P'] == 0.0 and r['L'] == 0.0:
            continue
        raw[r['epoch']][r['sat']][r['code']] = r

    # Determine obs types present per system
    sys_types = defaultdict(set)
    for ep_key, sat_dict in raw.items():
        for sat, code_dict in sat_dict.items():
            sys_ch = sat[0]
            for code2ch in code_dict:
                for ot in obs_code_to_rinex3(code2ch):
                    sys_types[sys_ch].add(ot)

    # Sort obs types by preferred order
    all_types = {}
    for sys_ch, types in sys_types.items():
        ordered = [t for t in OBS_ORDER if t in types]
        # Append any non-standard types not in OBS_ORDER
        ordered += sorted(t for t in types if t not in OBS_ORDER)
        all_types[sys_ch] = ordered

    # Build epoch list: {sat → {obs_type → (value, lli, snr_ind)}}
    epochs = []
    for ep_key in sorted(raw.keys()):
        sat_obs = {}
        for sat, code_dict in raw[ep_key].items():
            sys_ch = sat[0]
            obs_map = {}
            for code2ch, r in code_dict.items():
                snr_ind = snr_to_indicator(r['SNR'])
                lli     = r['LLI'] & 3
                if r['P'] != 0.0:
                    obs_map[f'C{code2ch}'] = (r['P'],   0,   0,   False)
                if r['L'] != 0.0:
                    obs_map[f'L{code2ch}'] = (r['L'],   lli, snr_ind, True)
                if r['D'] != 0.0:
                    obs_map[f'D{code2ch}'] = (r['D'],   0,   0,   False)
                if r['SNR'] > 0.0:
                    obs_map[f'S{code2ch}'] = (r['SNR'], 0,   0,   False)
            if obs_map:
                sat_obs[sat] = obs_map
        if sat_obs:
            epochs.append((ep_key, sat_obs))

    return epochs, all_types


# ── RINEX 3 writer ─────────────────────────────────────────────────────────────

def write_rinex3(path, epochs, all_types, pos_xyz, marker_name, log_basename):
    if not epochs:
        print('No observations to write.', file=sys.stderr)
        return

    first_ep = epochs[0][0]
    n_obs_total = sum(len(sd) for _, sd in epochs)

    # Collect system order (put GPS first)
    sys_order = sorted(all_types.keys(), key=lambda s: ('GRJECS').index(s)
                       if s in 'GRJECS' else 99)

    with open(path, 'w') as f:
        def w(line):
            f.write(f'{line:<60s}{line[60:] if len(line) > 60 else ""}\n'
                    if len(line) > 60
                    else f'{line:<60s}\n')

        def hdr(content, label):
            f.write(f'{content:<60s}{label:<20s}\n')

        # RINEX VERSION / TYPE
        hdr(f'{RINEX_VER:<9s}{"":11s}{"O":1s}{"":19s}{"M":1s}{"":19s}',
            'RINEX VERSION / TYPE')
        hdr(f'{PROG_NAME:<20s}{log_basename:<20s}', 'PGM / RUN BY / DATE')
        hdr(f'{marker_name:<60s}', 'MARKER NAME')
        hdr(f'{"ROVER":<20s}{"":40s}', 'MARKER TYPE')
        hdr(f'{"":20s}{"PocketSDR":20s}{"":20s}', 'REC # / TYPE / VERS')
        hdr(f'{"":20s}{"":20s}{"":20s}', 'ANT # / TYPE')

        # Approximate position
        if pos_xyz:
            x, y, z = pos_xyz
            hdr(f'{x:14.4f}{y:14.4f}{z:14.4f}{"":18s}', 'APPROX POSITION XYZ')
        else:
            hdr(f'{"0.0000":>14s}{"0.0000":>14s}{"0.0000":>14s}{"":18s}',
                'APPROX POSITION XYZ')

        hdr(f'{"1.0000":>14s}{"0.0000":>14s}{"0.0000":>14s}{"":18s}',
            'ANTENNA: DELTA H/E/N')

        # SYS / # / OBS TYPES (13 obs types per line max)
        for sys_ch in sys_order:
            types = all_types[sys_ch]
            n = len(types)
            first = True
            for start in range(0, n, 13):
                chunk = types[start:start+13]
                type_str = ''.join(f'{t:>4s}' for t in chunk)
                if first:
                    hdr(f'{sys_ch:1s}{n:>5d}{type_str:<52s}',
                        'SYS / # / OBS TYPES')
                    first = False
                else:
                    hdr(f'{"":6s}{type_str:<54s}', 'SYS / # / OBS TYPES')

        # TIME OF FIRST OBS
        y, mo, d, h, mi, s = first_ep
        si = int(s)
        frac = s - si
        hdr(f'{y:6d}{mo:6d}{d:6d}{h:6d}{mi:6d}{s:13.7f}{"":5s}{"GPS":3s}{"":9s}',
            'TIME OF FIRST OBS')

        hdr('', 'END OF HEADER')

        # ── data records ──────────────────────────────────────────────────────
        for ep_key, sat_dict in epochs:
            y, mo, d, h, mi, s = ep_key
            n_sat = len(sat_dict)
            f.write(f'> {y:4d} {mo:02d} {d:02d} {h:02d} {mi:02d}{s:11.7f}  0{n_sat:3d}\n')

            for sat in sorted(sat_dict.keys()):
                sys_ch = sat[0]
                if sys_ch not in all_types:
                    continue
                obs_types = all_types[sys_ch]
                obs_map   = sat_dict[sat]

                line = sat.ljust(3)
                for ot in obs_types:
                    if ot in obs_map:
                        val, lli, snr_ind, is_phase = obs_map[ot]
                        line += fmt_obs(val, lli, snr_ind, is_phase)
                    else:
                        line += ' ' * 16
                f.write(line.rstrip() + '\n')

    print(f'RINEX 3.03 OBS: {path}')
    print(f'  {len(epochs)} epochs, {n_obs_total} satellite-epochs, '
          f'systems: {" ".join(sys_order)}')
    for sys_ch in sys_order:
        print(f'  {sys_ch}: {" ".join(all_types[sys_ch])}')


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Convert pocket.log $OBS to RINEX 3.03 OBS for RTK post-processing')
    ap.add_argument('logfile', help='pocket.log input file')
    ap.add_argument('-o', '--out', default=None,
                    help='output RINEX OBS file (default: <logfile-stem>.obs)')
    ap.add_argument('--sys', default=None,
                    help='comma-separated satellite systems to include, e.g. G,J (default: all)')
    ap.add_argument('--marker', default='ROVER',
                    help='RINEX MARKER NAME (default: ROVER)')
    args = ap.parse_args()

    out_path = args.out or os.path.splitext(args.logfile)[0] + '.obs'
    sys_filter = set(args.sys.upper().split(',')) if args.sys else None

    print(f'Reading {args.logfile} ...', end='', flush=True)
    obs_records, pos_xyz = parse_log(args.logfile)
    print(f' {len(obs_records)} $OBS records')

    if not obs_records:
        print('No $OBS records found.', file=sys.stderr)
        sys.exit(1)

    epochs, all_types = group_epochs(obs_records, sys_filter)
    if not epochs:
        print('No valid observations after filtering.', file=sys.stderr)
        sys.exit(1)

    write_rinex3(out_path, epochs, all_types,
                 pos_xyz, args.marker, os.path.basename(args.logfile))

    print()
    if pos_xyz:
        print(f'Approx position from log: X={pos_xyz[0]:.1f} Y={pos_xyz[1]:.1f} Z={pos_xyz[2]:.1f}')
    else:
        print('No valid $POS found for approx position — header will have 0,0,0 (harmless).')

    print()
    print('Next steps:')
    print('  1. Download base CORS RINEX from Geoscience Australia:')
    print('       https://gnss.ga.gov.au/rinex/')
    print('     Nearest to Canberra: station CANB or TIDB')
    print('     Select date + hour window matching your experiment.')
    print()
    print('  2. Run RTK post-processing:')
    print(f'     rnx2rtkp -k <config> {out_path} base.obs base.nav -o rtk.pos')
    print()
    print('  3. Evaluate result:')
    print(f'     python3 python/pocket_pos_plot.py rtk.pos --ahd LAT LON AHD N --ant 0.5')


if __name__ == '__main__':
    main()
