#!/usr/bin/env python3
#
#  pocket_eph2rnx.py - Convert pocket.log $EPH lines to a RINEX 3.04 NAV file.
#
#  Broadcast ephemeris is identical for a given satellite/time window no
#  matter who decodes it -- it doesn't have to come from the base station.
#  pocket_trk already decodes and logs every ephemeris it receives ($EPH,
#  see src/sdr_pvt.c:out_log_eph), so the rover's own log is a complete,
#  immediately-available substitute for a base-station/CORS nav file --
#  no need to wait on GA (or any other provider) to publish theirs, which
#  can lag by a day or more after the UTC day ends.
#
#  Field layout mirrors src/sdr_pvt.c:out_log_eph() (GPS/QZS/GAL branch of
#  eph_t) and the writer mirrors lib/RTKLIB/src/rinex.c:outrnxnavb() exactly
#  (same D-exponent number format, same 8-line-per-record layout) so the
#  output is byte-format-compatible with what this repo's own RTKLIB would
#  produce.
#
#  Only GPS (G) and QZSS (J) are supported -- this covers the drone flight
#  config (script/real-time_ppp_l1l2_hybrid_drone.sh tracks GPS+QZS only).
#  GLONASS (R) uses a different record layout (geph_t/GNAV) and BeiDou (C)
#  needs a GPST->BDT conversion; both are skipped with a warning rather than
#  emitting untested/approximate records for constellations this flight
#  doesn't use.
#
#  Usage:
#    python3 pocket_eph2rnx.py pocket.log -o rover.nav
#
import sys, argparse, os, math
from collections import Counter
from datetime import datetime, timezone

GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc).timestamp()

URA_NOMINAL = [2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16.0, 32.0, 64.0, 128.0, 256.0,
               512.0, 1024.0, 2048.0, 4096.0, 8192.0]


def uravalue(sva):
    return URA_NOMINAL[sva] if 0 <= sva < 15 else 8192.0


def to_week_tow(ts):
    d = ts - GPS_EPOCH
    week = int(d // 604800)
    return week, d - week * 604800


def fmt_nav(value):
    """Mirror RTKLIB's outnavf(): ' 0.123456789012D+02', 12 mantissa digits."""
    e = 0 if abs(value) < 1e-99 else math.floor(math.log10(abs(value))) + 1
    mant = abs(value) / (10.0 ** (e - 12))
    sign = '-' if value < 0 else ' '
    return f' {sign}.{mant:012.0f}D{e:+03d}'


# ── parsing ────────────────────────────────────────────────────────────────────

EPH_FIELDS = ['iode', 'iodc', 'sva', 'svh', 'toe_ts', 'toc_ts', 'ttr_ts',
              'A', 'e', 'i0', 'OMG0', 'omg', 'M0', 'deln', 'OMGd', 'idot',
              'crc', 'crs', 'cuc', 'cus', 'cic', 'cis', 'toes', 'fit',
              'f0', 'f1', 'f2', 'tgd0', 'code', 'flag']
_INT_FIELDS = {'iode', 'iodc', 'sva', 'svh', 'toe_ts', 'toc_ts', 'ttr_ts',
               'code', 'flag'}


def parse_eph_log(path):
    records = {}  # (sat, toe_ts) -> field dict; keeps latest re-broadcast
    skipped = Counter()
    with open(path, errors='replace') as f:
        for line in f:
            if not line.startswith('$EPH,'):
                continue
            p = line.rstrip('\n').split(',')
            if len(p) < 34:
                continue
            sat = p[2]
            if sat[0] not in ('G', 'J'):
                skipped[sat[0]] += 1
                continue
            try:
                vals = [float(x) for x in p[4:34]]
            except ValueError:
                continue
            rec = dict(zip(EPH_FIELDS, vals))
            for k in _INT_FIELDS:
                rec[k] = int(rec[k])
            rec['sat'] = sat
            records[(sat, rec['toe_ts'])] = rec
    return sorted(records.values(), key=lambda r: (r['sat'], r['toe_ts'])), skipped


# ── RINEX 3.04 NAV writer (mirrors outrnxnavh/outrnxnavb) ──────────────────────

def write_rinex_nav(path, records, log_basename):
    with open(path, 'w') as f:
        f.write(f'{3.04:9.2f}{"":11s}{"N: GNSS NAV DATA":<20s}'
                f'{"M: Mixed":<20s}{"RINEX VERSION / TYPE":<20s}\n')
        date = datetime.now(timezone.utc).strftime('%Y%m%d %H%M%S UTC')
        f.write(f'{"pocket_eph2rnx.py":<20.20s}{log_basename:<20.20s}'
                f'{date:<20.20s}{"PGM / RUN BY / DATE":<20s}\n')
        f.write(f'{"":60s}{"END OF HEADER":<20s}\n')

        for r in records:
            sat = r['sat']
            sys = sat[0]
            ep = datetime.utcfromtimestamp(r['toc_ts'])
            f.write(f"{sat:<3s} {ep.year:04d} {ep.month:02d} {ep.day:02d} "
                    f"{ep.hour:02d} {ep.minute:02d} {ep.second:02d}")
            f.write(f"{fmt_nav(r['f0'])}{fmt_nav(r['f1'])}{fmt_nav(r['f2'])}\n")

            sep = '    '
            f.write(f"{sep}{fmt_nav(r['iode'])}{fmt_nav(r['crs'])}"
                    f"{fmt_nav(r['deln'])}{fmt_nav(r['M0'])}\n")
            f.write(f"{sep}{fmt_nav(r['cuc'])}{fmt_nav(r['e'])}"
                    f"{fmt_nav(r['cus'])}{fmt_nav(math.sqrt(r['A']))}\n")
            f.write(f"{sep}{fmt_nav(r['toes'])}{fmt_nav(r['cic'])}"
                    f"{fmt_nav(r['OMG0'])}{fmt_nav(r['cis'])}\n")
            f.write(f"{sep}{fmt_nav(r['i0'])}{fmt_nav(r['crc'])}"
                    f"{fmt_nav(r['omg'])}{fmt_nav(r['OMGd'])}\n")

            toe_week, _ = to_week_tow(r['toe_ts'])
            f.write(f"{sep}{fmt_nav(r['idot'])}{fmt_nav(r['code'])}"
                    f"{fmt_nav(toe_week)}{fmt_nav(r['flag'])}\n")

            f.write(f"{sep}{fmt_nav(uravalue(r['sva']))}{fmt_nav(r['svh'])}"
                    f"{fmt_nav(r['tgd0'])}{fmt_nav(r['iodc'])}\n")

            ttr_week, ttr_tow = to_week_tow(r['ttr_ts'])
            ttr_adj = ttr_tow + (ttr_week - toe_week) * 604800.0
            fit = 1.0 if (sys == 'J' and r['fit'] > 2.0) else \
                (0.0 if sys == 'J' else r['fit'])
            f.write(f"{sep}{fmt_nav(ttr_adj)}{fmt_nav(fit)}\n")

    print(f'RINEX 3.04 NAV: {path}  ({len(records)} ephemerides)')
    by_sys = Counter(r['sat'][0] for r in records)
    for s, n in sorted(by_sys.items()):
        print(f'  {s}: {n}')


def main():
    ap = argparse.ArgumentParser(
        description='Convert pocket.log $EPH to a RINEX 3.04 NAV file '
                    '(GPS+QZS only)')
    ap.add_argument('logfile', help='pocket.log input file')
    ap.add_argument('-o', '--out', default=None,
                    help='output RINEX NAV file (default: <logfile-stem>.nav)')
    args = ap.parse_args()

    out_path = args.out or os.path.splitext(args.logfile)[0] + '.nav'

    print(f'Reading {args.logfile} ...', end='', flush=True)
    records, skipped = parse_eph_log(args.logfile)
    print(f' {len(records)} unique ephemerides (GPS+QZS)')

    if skipped:
        skip_str = ', '.join(f'{s}:{n}' for s, n in sorted(skipped.items()))
        print(f'Skipped (unsupported system): {skip_str}', file=sys.stderr)

    if not records:
        print('No usable $EPH records found -- was the log written with '
              '-log and trace level >= 3?', file=sys.stderr)
        sys.exit(1)

    write_rinex_nav(out_path, records, os.path.basename(args.logfile))


if __name__ == '__main__':
    main()
