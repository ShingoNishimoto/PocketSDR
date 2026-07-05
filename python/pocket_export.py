#!/usr/bin/env python3
#
#  pocket_export.py - Export PocketSDR pocket.log to geo/data formats.
#
#  Reads $POS (and optionally $REFPOS) entries and writes any combination of:
#    KML     (Google Earth, PPP and SPP as separate folders)
#    GPX     (GPS track format, PPP and SPP as separate tracks)
#    GeoJSON (web mapping, per-point properties)
#    CSV     (spreadsheet-compatible)
#    Excel   (.xlsx with All / PPP / SPP sheets; Ref / RefPPP / RefSPP added when
#             $REFPOS records are present)
#
#  Usage:
#    python3 pocket_export.py pocket.log [--all] [--kml F] [--gpx F]
#                             [--geojson F] [--csv F] [--excel F]
#                             [--ref-kml F] [--ref-gpx F] [--ref-geojson F]
#                             [--ref-csv F] [--ref-excel F]
#
#  With --all, all formats are written; file names are derived from the log
#  file base name (e.g. pocket.log → pocket.kml, pocket.gpx, ...).
#  Reference solution files are written automatically when $REFPOS records
#  exist (pocket_ref.kml, pocket_ref.csv, …).
#
import sys, argparse, json, os
from datetime import datetime, timezone, timedelta
import numpy as np
import pandas as pd

# WGS84 ellipsoid parameters
_WGS84_A  = 6378137.0
_WGS84_E2 = 2 / 298.257223563 - (1 / 298.257223563) ** 2

SOLQ_LABEL = {0:'NONE', 1:'FIX', 2:'FLOAT', 3:'SBAS', 4:'DGPS',
              5:'SPP', 6:'PPP', 7:'DR'}

# GPS is ahead of UTC by 18 leap seconds as of 2017-01-01; update if needed.
GPS_LEAP_SECONDS = 18

# ── parsing ────────────────────────────────────────────────────────────────────

def parse_log(path):
    """Return list of dicts parsed from $POS lines in a pocket.log file.

    $POS format:
      t,year,month,day,hour,min,sec,lat,lon,hgt,q,ns,stdn,stde,stdu,dtr [17 fields]
      + ecef_x,ecef_y,ecef_z,vel_x,vel_y,vel_z                          [+6, len>=23]
      + gdop,pdop,hdop,vdop                                              [+4, len>=27]
      + dtr_drift                                                        [+1, len>=28]
    """
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$POS,'):
                continue
            p = line.split(',')
            if len(p) < 17:
                continue
            try:
                r = {
                    't':     float(p[1]),
                    'year':  int(float(p[2])), 'month': int(float(p[3])),
                    'day':   int(float(p[4])), 'hour':  int(float(p[5])),
                    'min':   int(float(p[6])), 'sec':   float(p[7]),
                    'lat':   float(p[8]),  'lon':  float(p[9]),
                    'hgt':   float(p[10]),
                    'q':     int(p[11]),   'ns':   int(p[12]),
                    'stdn':  float(p[13]), 'stde': float(p[14]),
                    'stdu':  float(p[15]), 'dtr':  float(p[16]),
                }
                if len(p) >= 23:
                    r['ecef_x'] = float(p[17])
                    r['ecef_y'] = float(p[18])
                    r['ecef_z'] = float(p[19])
                    r['vel_x']  = float(p[20])
                    r['vel_y']  = float(p[21])
                    r['vel_z']  = float(p[22])
                else:
                    r['ecef_x'], r['ecef_y'], r['ecef_z'] = \
                        llh_to_ecef(r['lat'], r['lon'], r['hgt'])
                    r['vel_x'] = r['vel_y'] = r['vel_z'] = float('nan')
                if len(p) >= 27:
                    r['gdop'] = float(p[23])
                    r['pdop'] = float(p[24])
                    r['hdop'] = float(p[25])
                    r['vdop'] = float(p[26])
                if len(p) >= 28:
                    r['dtr_drift'] = float(p[27])
                rows.append(r)
            except (ValueError, IndexError):
                continue
    return rows


def parse_sc_aowr_log(path):
    """Parse $LOG,<t>,SC_AOWR entries and return a dict keyed by t.

    Line format:  $LOG,<t>,SC_AOWR clk=<float> drift=<float> dt=<float> hist=<int>
    Returns: {t: {'clk': float, 'drift': float, 'dt': float, 'hist': int}, ...}
    """
    result = {}
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$LOG,'):
                continue
            p = line.split(',', 2)
            if len(p) < 3 or not p[2].startswith('SC_AOWR '):
                continue
            try:
                t     = float(p[1])
                msg   = p[2]
                clk   = float(msg.split('clk=')[1].split()[0])
                drift = float(msg.split('drift=')[1].split()[0])
                dt    = float(msg.split('dt=')[1].split()[0])
                hist  = int(msg.split('hist=')[1].split()[0])
                result[t] = {'clk': clk, 'drift': drift, 'dt': dt, 'hist': hist}
            except (ValueError, IndexError):
                continue
    return result


def parse_pppos_dbg_log(path):
    """Parse $LOG,<t>,PPPOS_DBG entries for PPP convergence monitoring.

    Handles both old format (no clk_std/zwd/amb_std fields) and new format.
    Returns: {t: {pos_std, stat, ndual, clk, clk_std, zwd, zwd_std,
                  amb_std, nconv, res_c}, ...}
    """
    def _fv(msg, key):
        try:
            return float(msg.split(key + '=')[1].split()[0])
        except (IndexError, ValueError):
            return float('nan')

    def _iv(msg, key):
        try:
            return int(float(msg.split(key + '=')[1].split()[0]))
        except (IndexError, ValueError):
            return None

    result = {}
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$LOG,'):
                continue
            p = line.split(',', 2)
            if len(p) < 3 or not p[2].startswith('PPPOS_DBG '):
                continue
            try:
                t   = float(p[1])
                msg = p[2]
                result[t] = {
                    'pos_std': _fv(msg, 'pos_std'),
                    'stat':    _iv(msg, 'stat'),
                    'ndual':   _iv(msg, 'ndual'),
                    'clk':     _fv(msg, 'clk'),
                    'clk_std': _fv(msg, 'clk_std'),
                    'zwd':     _fv(msg, 'zwd'),
                    'zwd_std': _fv(msg, 'zwd_std'),
                    'amb_std': _fv(msg, 'amb_std'),
                    'nconv':   _iv(msg, 'nconv'),
                    'res_c':   _fv(msg, 'res_c'),
                }
            except (ValueError, IndexError):
                continue
    return result


def parse_spp_seed_log(path, tag):
    """Parse $LOG entries with a given tag and return a dict keyed by t.

    Handles tags emitted by sdr_pvt.c:
      SPP_SEED     – main SPP seed (AOWR-corrected obs; dtr=0 when clock fixed)
      REF_SPP_SEED – reference SPP seed in PPP mode (uncorrected obs, free clock)
      REF_SPP      – reference SPP in SPP mode (uncorrected obs, free clock)

    Line format:  $LOG,<t>,<tag> stat=<int> dtr=<float> [msg=<str>]
    Returns: {t: {'stat': int, 'dtr': float}, ...}
    """
    result = {}
    prefix = tag + ' '
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$LOG,'):
                continue
            p = line.split(',', 2)
            if len(p) < 3 or not p[2].startswith(prefix):
                continue
            try:
                t    = float(p[1])
                msg  = p[2]
                stat = int(msg.split('stat=')[1].split()[0])
                dtr  = float(msg.split('dtr=')[1].split()[0])
                result[t] = {'stat': stat, 'dtr': dtr}
            except (ValueError, IndexError):
                continue
    return result


def parse_refpos_log(path):
    """Return list of dicts parsed from $REFPOS lines (reference/uncorrected solver).

    $REFPOS format (4-unknown solve, clock estimated):
      t,year,month,day,hour,min,sec,lat,lon,hgt,stat,ns,stdn,stde,stdu,dtr [17 fields]
      + gdop,pdop,hdop,vdop                                                 [+4, len>=21]
      + dtr_drift                                                           [+1, len>=22]
      + vel_x,vel_y,vel_z (ECEF m/s)                                       [+3, len>=25]
    ECEF position is computed from lat/lon/hgt (not in the log format).
    """
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$REFPOS,'):
                continue
            p = line.split(',')
            if len(p) < 17:
                continue
            try:
                r = {
                    't':     float(p[1]),
                    'year':  int(float(p[2])), 'month': int(float(p[3])),
                    'day':   int(float(p[4])), 'hour':  int(float(p[5])),
                    'min':   int(float(p[6])), 'sec':   float(p[7]),
                    'lat':   float(p[8]),  'lon':  float(p[9]),
                    'hgt':   float(p[10]),
                    'q':     int(p[11]),   'ns':   int(p[12]),
                    'stdn':  float(p[13]), 'stde': float(p[14]),
                    'stdu':  float(p[15]), 'dtr':  float(p[16]),
                }
                r['ecef_x'], r['ecef_y'], r['ecef_z'] = \
                    llh_to_ecef(r['lat'], r['lon'], r['hgt'])
                r['vel_x'] = r['vel_y'] = r['vel_z'] = float('nan')
                if len(p) >= 21:
                    r['gdop'] = float(p[17])
                    r['pdop'] = float(p[18])
                    r['hdop'] = float(p[19])
                    r['vdop'] = float(p[20])
                if len(p) >= 22:
                    r['dtr_drift'] = float(p[21])
                if len(p) >= 25:
                    r['vel_x'] = float(p[22])
                    r['vel_y'] = float(p[23])
                    r['vel_z'] = float(p[24])
                rows.append(r)
            except (ValueError, IndexError):
                continue
    return rows


def gpst_to_utc(r):
    """Convert $POS/$REFPOS epoch fields (GPST) to UTC datetime object."""
    si = int(r['sec'])
    us = int((r['sec'] - si) * 1e6)
    # time2epoch() can return sec=60.000 at a minute boundary due to float rounding
    extra_min, si = divmod(si, 60)
    dt = datetime(r['year'], r['month'], r['day'],
                  r['hour'], r['min'], si, us, tzinfo=timezone.utc)
    return dt + timedelta(minutes=extra_min) - timedelta(seconds=GPS_LEAP_SECONDS)

def gpst_str(r):
    dt = gpst_to_utc(r) + timedelta(seconds=GPS_LEAP_SECONDS)
    frac = r['sec'] % 1
    return (f"{dt.year}-{dt.month:02d}-{dt.day:02d}"
            f"T{dt.hour:02d}:{dt.minute:02d}:{dt.second + frac:06.3f}")

def llh_to_ecef(lat_deg, lon_deg, hgt_m):
    """WGS84 geodetic (deg, deg, m) → ECEF (m)."""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    N = _WGS84_A / np.sqrt(1 - _WGS84_E2 * np.sin(lat) ** 2)
    X = (N + hgt_m) * np.cos(lat) * np.cos(lon)
    Y = (N + hgt_m) * np.cos(lat) * np.sin(lon)
    Z = (N * (1 - _WGS84_E2) + hgt_m) * np.sin(lat)
    return X, Y, Z


# ── KML ───────────────────────────────────────────────────────────────────────

def write_kml(rows, path, label=''):
    ppp = [r for r in rows if r['q'] == 6]
    spp = [r for r in rows if r['q'] != 6]
    tag = f' ({label})' if label else ''

    def track(pts, name, line_style, dot_style):
        coords = ' '.join(f"{r['lon']},{r['lat']},{r['hgt']}" for r in pts)
        pms = []
        if coords:
            pms.append(
                f'    <Placemark><name>{name} Track</name>'
                f'<styleUrl>#{line_style}</styleUrl>'
                f'<LineString><altitudeMode>absolute</altitudeMode>'
                f'<coordinates>{coords}</coordinates>'
                f'</LineString></Placemark>')
        for r in pts:
            desc = (f"t={r['t']:.1f}s Q={SOLQ_LABEL.get(r['q'],r['q'])} "
                    f"ns={r['ns']} σN={r['stdn']:.3f}m σE={r['stde']:.3f}m "
                    f"σU={r['stdu']:.3f}m")
            pms.append(
                f'    <Placemark><description>{desc}</description>'
                f'<styleUrl>#{dot_style}</styleUrl>'
                f'<Point><altitudeMode>absolute</altitudeMode>'
                f'<coordinates>{r["lon"]},{r["lat"]},{r["hgt"]}</coordinates>'
                f'</Point></Placemark>')
        return '\n'.join(pms)

    kml = f'''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>PocketSDR Position{tag}</name>
  <Style id="pppLine"><LineStyle><color>ff0000ff</color><width>2</width></LineStyle></Style>
  <Style id="pppDot"><IconStyle><color>ff0000ff</color><scale>0.5</scale></IconStyle>
                    <LabelStyle><scale>0</scale></LabelStyle></Style>
  <Style id="sppLine"><LineStyle><color>ff00ff00</color><width>1</width></LineStyle></Style>
  <Style id="sppDot"><IconStyle><color>ff00ff00</color><scale>0.4</scale></IconStyle>
                    <LabelStyle><scale>0</scale></LabelStyle></Style>
  <Folder>
    <name>PPP (Q=6) — {len(ppp)} epochs</name>
{track(ppp, "PPP", "pppLine", "pppDot")}
  </Folder>
  <Folder>
    <name>SPP (Q=5) — {len(spp)} epochs</name>
{track(spp, "SPP", "sppLine", "sppDot")}
  </Folder>
</Document>
</kml>'''
    with open(path, 'w') as f:
        f.write(kml)
    kind = f'KML{tag}:'
    print(f'{kind:<16} {path}  ({len(ppp)} PPP + {len(spp)} SPP)')

# ── GPX ───────────────────────────────────────────────────────────────────────

def write_gpx(rows, path, label=''):
    def trkpt(r):
        dt = gpst_to_utc(r)
        ts = dt.strftime('%Y-%m-%dT%H:%M:%S.') + f'{dt.microsecond//1000:03d}Z'
        ql = SOLQ_LABEL.get(r['q'], str(r['q']))
        return (f'      <trkpt lat="{r["lat"]:.9f}" lon="{r["lon"]:.9f}">'
                f'<ele>{r["hgt"]:.3f}</ele><time>{ts}</time>'
                f'<desc>{ql} ns={r["ns"]} σN={r["stdn"]:.3f}m</desc>'
                f'</trkpt>')

    ppp_pts = '\n'.join(trkpt(r) for r in rows if r['q'] == 6)
    spp_pts = '\n'.join(trkpt(r) for r in rows if r['q'] != 6)
    ppp_n = sum(1 for r in rows if r['q'] == 6)
    spp_n = len(rows) - ppp_n
    tag = f' ({label})' if label else ''

    gpx = f'''<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="PocketSDR pocket_export.py"
     xmlns="http://www.topografix.com/GPX/1/1">
  <trk>
    <name>PPP (Q=6){tag}</name>
    <trkseg>
{ppp_pts}
    </trkseg>
  </trk>
  <trk>
    <name>SPP (Q=5){tag}</name>
    <trkseg>
{spp_pts}
    </trkseg>
  </trk>
</gpx>'''
    with open(path, 'w') as f:
        f.write(gpx)
    kind = f'GPX{tag}:'
    print(f'{kind:<16} {path}  ({ppp_n} PPP + {spp_n} SPP)')

# ── GeoJSON ───────────────────────────────────────────────────────────────────

def write_geojson(rows, path, label=''):
    source = label or 'main'
    features = [{
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': [round(r['lon'], 9), round(r['lat'], 9),
                            round(r['hgt'], 3)],
        },
        'properties': {
            'source':         source,
            't_s':            r['t'],
            'datetime_gpst':  gpst_str(r),
            'quality':        r['q'],
            'quality_label':  SOLQ_LABEL.get(r['q'], str(r['q'])),
            'n_sats':         r['ns'],
            'stdn_m':         round(r['stdn'], 4),
            'stde_m':         round(r['stde'], 4),
            'stdu_m':         round(r['stdu'], 4),
            'dtr_s':          round(r['dtr'], 9),
        },
    } for r in rows]
    fc = {'type': 'FeatureCollection', 'features': features}
    with open(path, 'w') as f:
        json.dump(fc, f, separators=(',', ':'))
    tag = f' ({label})' if label else ''
    kind = f'GeoJSON{tag}:'
    print(f'{kind:<16} {path}  ({len(features)} features)')

# ── CSV / Excel ───────────────────────────────────────────────────────────────

def make_df(rows):
    """Build a DataFrame from $POS or $REFPOS row dicts.

    Optional fields (ecef, velocity, DOP, SPP seed clocks) use .get() so the
    same function works for both row types; missing values become NaN.
    SPP seed clock fields are merged into row dicts by main() before this call:
      spp_seed_stat / spp_seed_dtr_s  – main SPP seed (AOWR-corrected)
      ref_spp_stat  / ref_spp_dtr_s   – reference SPP seed (uncorrected, free clock)
    """
    _nan = float('nan')
    return pd.DataFrame([{
        't_s':           r['t'],
        'datetime_gpst': gpst_str(r),
        # Geodetic
        'lat_deg':       r['lat'],
        'lon_deg':       r['lon'],
        'hgt_m':         r['hgt'],
        # ECEF position (NaN for $REFPOS rows which don't carry ECEF)
        'ecef_x_m':      r.get('ecef_x', _nan),
        'ecef_y_m':      r.get('ecef_y', _nan),
        'ecef_z_m':      r.get('ecef_z', _nan),
        # ECEF velocity (non-zero only in kinematic PPP; NaN for $REFPOS)
        'vel_x_m_s':     r.get('vel_x', _nan),
        'vel_y_m_s':     r.get('vel_y', _nan),
        'vel_z_m_s':     r.get('vel_z', _nan),
        # Solution quality
        'quality':       r['q'],
        'quality_label': SOLQ_LABEL.get(r['q'], str(r['q'])),
        'n_sats':        r['ns'],
        'stdn_m':        r['stdn'],
        'stde_m':        r['stde'],
        'stdu_m':        r['stdu'],
        'dtr_s':         r['dtr'],
        # Clock drift rate (s/s): Doppler-derived in SPP; WLS in SC AOWR; 0 in PPP
        'dtr_drift_s':   r.get('dtr_drift', _nan),
        # DOP (NaN for old logs without DOP fields)
        # For $POS (clock-corrected): gdop=0, pdop=PDOP3 (position-only, no clock coupling)
        # For $REFPOS (4-unknown):    gdop=GDOP, pdop=PDOP with clock-position coupling
        'gdop':          r.get('gdop', _nan),
        'pdop':          r.get('pdop', _nan),
        'hdop':          r.get('hdop', _nan),
        'vdop':          r.get('vdop', _nan),
        # SPP seed receiver clock offsets (debug; NaN when not in -ps_sc mode or no fix)
        'spp_seed_stat':  r.get('spp_seed_stat', _nan),
        'spp_seed_dtr_s': r.get('spp_seed_dtr',  _nan),
        'ref_spp_stat':   r.get('ref_spp_stat',  _nan),
        'ref_spp_dtr_s':  r.get('ref_spp_dtr',   _nan),
        # SC_AOWR clock estimate and drift rate (NaN when AOWR not yet active)
        'sc_aowr_clk_s':  r.get('sc_aowr_clk',   _nan),
        'sc_aowr_drift':  r.get('sc_aowr_drift',  _nan),
        # PPP KF convergence diagnostics from PPPOS_DBG (NaN when not in PPP mode
        # or old logs without these fields)
        'ppp_kf_pos_std_m': r.get('ppp_pos_std', _nan),
        'ppp_kf_ndual':     r.get('ppp_ndual',   _nan),
        'ppp_clk_m':        r.get('ppp_clk',     _nan),
        'ppp_clk_std_m':    r.get('ppp_clk_std', _nan),
        'ppp_zwd_m':        r.get('ppp_zwd',     _nan),
        'ppp_zwd_std_m':    r.get('ppp_zwd_std', _nan),
        'ppp_amb_std_m':    r.get('ppp_amb_std', _nan),
        'ppp_nconv':        r.get('ppp_nconv',   _nan),
        'ppp_res_c_m':      r.get('ppp_res_c',   _nan),
    } for r in rows])

def write_csv(rows, path, label=''):
    df = make_df(rows)
    df.to_csv(path, index=False, float_format='%.9g')
    tag = f' ({label})' if label else ''
    kind = f'CSV{tag}:'
    print(f'{kind:<16} {path}  ({len(df)} rows)')

def write_excel(rows, path, ref_rows=None):
    df = make_df(rows)
    if 'quality' in df.columns:
        ppp_df = df[df['quality'] == 6]
        spp_df = df[df['quality'] != 6]
    else:
        ppp_df = spp_df = df
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        df.to_excel(xw, sheet_name='All', index=False)
        ppp_df.to_excel(xw, sheet_name='PPP', index=False)
        spp_df.to_excel(xw, sheet_name='SPP', index=False)
        if ref_rows:
            rdf = make_df(ref_rows)
            rppp = rdf[rdf['quality'] == 6]
            rspp = rdf[rdf['quality'] != 6]
            rdf.to_excel(xw, sheet_name='Ref', index=False)
            rppp.to_excel(xw, sheet_name='RefPPP', index=False)
            rspp.to_excel(xw, sheet_name='RefSPP', index=False)
    ref_note = f' + {len(ref_rows)} ref' if ref_rows else ''
    print(f'Excel:           {path}  ({len(df)} rows{ref_note})')

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Export PocketSDR log to KML / GPX / GeoJSON / CSV / Excel.\n'
                    '$REFPOS records (SC AOWR reference solver) are exported '
                    'automatically with --all, or via --ref-* flags.')
    ap.add_argument('logfile', help='pocket.log file')
    ap.add_argument('--all', action='store_true',
                    help='write all formats; names derived from logfile stem')
    # Main ($POS) output
    ap.add_argument('--kml',     metavar='FILE')
    ap.add_argument('--gpx',     metavar='FILE')
    ap.add_argument('--geojson', metavar='FILE')
    ap.add_argument('--csv',     metavar='FILE')
    ap.add_argument('--excel',   metavar='FILE')
    # Reference ($REFPOS) output
    ap.add_argument('--ref-kml',     metavar='FILE')
    ap.add_argument('--ref-gpx',     metavar='FILE')
    ap.add_argument('--ref-geojson', metavar='FILE')
    ap.add_argument('--ref-csv',     metavar='FILE')
    args = ap.parse_args()

    rows = parse_log(args.logfile)
    if rows:
        ppp_n = sum(1 for r in rows if r['q'] == 6)
        spp_n = len(rows) - ppp_n
        print(f'Parsed {len(rows)} $POS records from {args.logfile}'
              f'  ({ppp_n} PPP, {spp_n} SPP)')
    else:
        print(f'No $POS records in {args.logfile}', file=sys.stderr)

    ref_rows = parse_refpos_log(args.logfile)
    if ref_rows:
        rppp_n = sum(1 for r in ref_rows if r['q'] == 6)
        rspp_n = len(ref_rows) - rppp_n
        print(f'Parsed {len(ref_rows)} $REFPOS records'
              f'  ({rppp_n} PPP, {rspp_n} SPP)')

    if not rows and not ref_rows:
        print('No position records found, exiting.', file=sys.stderr)
        sys.exit(1)

    # Merge SPP seed clock entries into rows for debug export.
    # REF_SPP_SEED (PPP mode) and REF_SPP (SPP mode) both carry the reference
    # solver clock; combine them into a single lookup.
    spp_clk = parse_spp_seed_log(args.logfile, 'SPP_SEED')
    ref_spp_clk = {**parse_spp_seed_log(args.logfile, 'REF_SPP_SEED'),
                   **parse_spp_seed_log(args.logfile, 'REF_SPP')}
    if spp_clk or ref_spp_clk:
        for r in rows:
            if r['t'] in spp_clk:
                r['spp_seed_stat'] = spp_clk[r['t']]['stat']
                r['spp_seed_dtr']  = spp_clk[r['t']]['dtr']
            if r['t'] in ref_spp_clk:
                r['ref_spp_stat'] = ref_spp_clk[r['t']]['stat']
                r['ref_spp_dtr']  = ref_spp_clk[r['t']]['dtr']
        print(f'Merged {len(spp_clk)} SPP_SEED and {len(ref_spp_clk)} REF_SPP clock entries')

    sc_aowr = parse_sc_aowr_log(args.logfile)
    if sc_aowr:
        for r in rows:
            if r['t'] in sc_aowr:
                r['sc_aowr_clk']   = sc_aowr[r['t']]['clk']
                r['sc_aowr_drift'] = sc_aowr[r['t']]['drift']
        print(f'Merged {len(sc_aowr)} SC_AOWR clock/drift entries')

    pppos_dbg = parse_pppos_dbg_log(args.logfile)
    if pppos_dbg:
        for r in rows:
            if r['t'] in pppos_dbg:
                d = pppos_dbg[r['t']]
                r['ppp_pos_std'] = d['pos_std']
                r['ppp_ndual']   = d['ndual']
                r['ppp_clk']     = d['clk']
                r['ppp_clk_std'] = d['clk_std']
                r['ppp_zwd']     = d['zwd']
                r['ppp_zwd_std'] = d['zwd_std']
                r['ppp_amb_std'] = d['amb_std']
                r['ppp_nconv']   = d['nconv']
                r['ppp_res_c']   = d['res_c']
        print(f'Merged {len(pppos_dbg)} PPPOS_DBG convergence entries')

    stem = os.path.splitext(args.logfile)[0]
    ref_stem = stem + '_ref'

    want = args.all or any([args.kml, args.gpx, args.geojson, args.csv, args.excel,
                            args.ref_kml, args.ref_gpx, args.ref_geojson, args.ref_csv])
    if not want:
        ap.print_help()
        sys.exit(0)

    # ── main ($POS) output ──
    if args.all or args.kml:
        write_kml(rows,     args.kml     or stem + '.kml')
    if args.all or args.gpx:
        write_gpx(rows,     args.gpx     or stem + '.gpx')
    if args.all or args.geojson:
        write_geojson(rows, args.geojson or stem + '.geojson')
    if args.all or args.csv:
        write_csv(rows,     args.csv     or stem + '.csv')
    if args.all or args.excel:
        # Include ref sheets in the same Excel file for easy comparison
        write_excel(rows,   args.excel   or stem + '.xlsx',
                    ref_rows=ref_rows or None)

    # ── reference ($REFPOS) output ──
    if ref_rows:
        if args.all or args.ref_kml:
            write_kml(ref_rows,     args.ref_kml     or ref_stem + '.kml',
                      label='ref')
        if args.all or args.ref_gpx:
            write_gpx(ref_rows,     args.ref_gpx     or ref_stem + '.gpx',
                      label='ref')
        if args.all or args.ref_geojson:
            write_geojson(ref_rows, args.ref_geojson or ref_stem + '.geojson',
                          label='ref')
        if args.all or args.ref_csv:
            write_csv(ref_rows,     args.ref_csv     or ref_stem + '.csv',
                      label='ref')
    elif any([args.ref_kml, args.ref_gpx, args.ref_geojson, args.ref_csv]):
        print('Warning: no $REFPOS records found — ref output skipped',
              file=sys.stderr)

if __name__ == '__main__':
    main()
