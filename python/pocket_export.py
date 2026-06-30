#!/usr/bin/env python3
#
#  pocket_export.py - Export PocketSDR pocket.log to geo/data formats.
#
#  Reads $POS entries and writes any combination of:
#    KML     (-  Google Earth, PPP and SPP as separate folders)
#    GPX     (- GPS track format, PPP and SPP as separate tracks)
#    GeoJSON (- web mapping, per-point properties)
#    CSV     (- spreadsheet-compatible)
#    Excel   (- .xlsx with All / PPP / SPP sheets)
#
#  Usage:
#    python3 pocket_export.py pocket.log [--all] [--kml F] [--gpx F]
#                             [--geojson F] [--csv F] [--excel F]
#
#  With --all, all formats are written; file names are derived from the log
#  file base name (e.g. pocket.log → pocket.kml, pocket.gpx, ...).
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

    $POS format (fields 17-22 added in newer builds):
      t,year,month,day,hour,min,sec,lat,lon,hgt,q,ns,stdn,stde,stdu,dtr,
      ecef_x,ecef_y,ecef_z,vel_x,vel_y,vel_z
    Older logs with only 17 fields fall back to computed ECEF and zero velocity.
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
                    # New format: ECEF position and velocity from RTKLIB
                    r['ecef_x'] = float(p[17])
                    r['ecef_y'] = float(p[18])
                    r['ecef_z'] = float(p[19])
                    r['vel_x']  = float(p[20])
                    r['vel_y']  = float(p[21])
                    r['vel_z']  = float(p[22])
                else:
                    # Old format: derive ECEF from lat/lon/hgt, velocity unknown
                    r['ecef_x'], r['ecef_y'], r['ecef_z'] = \
                        llh_to_ecef(r['lat'], r['lon'], r['hgt'])
                    r['vel_x'] = r['vel_y'] = r['vel_z'] = float('nan')
                rows.append(r)
            except (ValueError, IndexError):
                continue
    return rows

def gpst_to_utc(r):
    """Convert $POS epoch fields (GPST) to UTC datetime object."""
    si = int(r['sec'])
    us = int((r['sec'] - si) * 1e6)
    # time2epoch() can return sec=60.000 at a minute boundary due to float rounding
    extra_min, si = divmod(si, 60)
    dt = datetime(r['year'], r['month'], r['day'],
                  r['hour'], r['min'], si, us, tzinfo=timezone.utc)
    return dt + timedelta(minutes=extra_min) - timedelta(seconds=GPS_LEAP_SECONDS)

def gpst_str(r):
    # Use gpst_to_utc then add back leap seconds — handles sec=60 rollover cleanly
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

def write_kml(rows, path):
    ppp = [r for r in rows if r['q'] == 6]
    spp = [r for r in rows if r['q'] != 6]

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
  <name>PocketSDR Position</name>
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
    print(f'KML:     {path}  ({len(ppp)} PPP + {len(spp)} SPP)')

# ── GPX ───────────────────────────────────────────────────────────────────────

def write_gpx(rows, path):
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

    gpx = f'''<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="PocketSDR pocket_export.py"
     xmlns="http://www.topografix.com/GPX/1/1">
  <trk>
    <name>PPP (Q=6)</name>
    <trkseg>
{ppp_pts}
    </trkseg>
  </trk>
  <trk>
    <name>SPP (Q=5)</name>
    <trkseg>
{spp_pts}
    </trkseg>
  </trk>
</gpx>'''
    with open(path, 'w') as f:
        f.write(gpx)
    print(f'GPX:     {path}  ({ppp_n} PPP + {spp_n} SPP)')

# ── GeoJSON ───────────────────────────────────────────────────────────────────

def write_geojson(rows, path):
    features = [{
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': [round(r['lon'], 9), round(r['lat'], 9),
                            round(r['hgt'], 3)],
        },
        'properties': {
            't_s': r['t'],
            'datetime_gpst': gpst_str(r),
            'quality': r['q'],
            'quality_label': SOLQ_LABEL.get(r['q'], str(r['q'])),
            'n_sats': r['ns'],
            'stdn_m': round(r['stdn'], 4),
            'stde_m': round(r['stde'], 4),
            'stdu_m': round(r['stdu'], 4),
            'dtr_s':  round(r['dtr'], 9),
        },
    } for r in rows]
    fc = {'type': 'FeatureCollection', 'features': features}
    with open(path, 'w') as f:
        json.dump(fc, f, separators=(',', ':'))
    print(f'GeoJSON: {path}  ({len(features)} features)')

# ── CSV / Excel ───────────────────────────────────────────────────────────────

def make_df(rows):
    return pd.DataFrame([{
        't_s':           r['t'],
        'datetime_gpst': gpst_str(r),
        # Geodetic
        'lat_deg':       r['lat'],
        'lon_deg':       r['lon'],
        'hgt_m':         r['hgt'],
        # ECEF position (converted from lat/lon/hgt via WGS84)
        'ecef_x_m':      r['ecef_x'],
        'ecef_y_m':      r['ecef_y'],
        'ecef_z_m':      r['ecef_z'],
        # ECEF velocity (numerical derivative; meaningful only for kinematic mode)
        'vel_x_m_s':     r['vel_x'],
        'vel_y_m_s':     r['vel_y'],
        'vel_z_m_s':     r['vel_z'],
        # Solution quality
        'quality':       r['q'],
        'quality_label': SOLQ_LABEL.get(r['q'], str(r['q'])),
        'n_sats':        r['ns'],
        'stdn_m':        r['stdn'],
        'stde_m':        r['stde'],
        'stdu_m':        r['stdu'],
        'dtr_s':         r['dtr'],
    } for r in rows])

def write_csv(rows, path):
    df = make_df(rows)
    df.to_csv(path, index=False, float_format='%.9g')
    print(f'CSV:     {path}  ({len(df)} rows)')

def write_excel(rows, path):
    df = make_df(rows)
    ppp_df = df[df['quality'] == 6]
    spp_df = df[df['quality'] != 6]
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        df.to_excel(xw, sheet_name='All', index=False)
        ppp_df.to_excel(xw, sheet_name='PPP', index=False)
        spp_df.to_excel(xw, sheet_name='SPP', index=False)
    print(f'Excel:   {path}  ({len(df)} rows: {len(ppp_df)} PPP, {len(spp_df)} SPP)')

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Export PocketSDR log to KML / GPX / GeoJSON / CSV / Excel')
    ap.add_argument('logfile', help='pocket.log file')
    ap.add_argument('--all', action='store_true',
                    help='write all formats; names derived from logfile stem')
    ap.add_argument('--kml',     metavar='FILE')
    ap.add_argument('--gpx',     metavar='FILE')
    ap.add_argument('--geojson', metavar='FILE')
    ap.add_argument('--csv',     metavar='FILE')
    ap.add_argument('--excel',   metavar='FILE')
    args = ap.parse_args()

    rows = parse_log(args.logfile)
    if not rows:
        print(f'No $POS records in {args.logfile}', file=sys.stderr)
        sys.exit(1)

    ppp_n = sum(1 for r in rows if r['q'] == 6)
    spp_n = len(rows) - ppp_n
    print(f'Parsed {len(rows)} $POS records from {args.logfile}'
          f'  ({ppp_n} PPP, {spp_n} SPP)')

    stem = os.path.splitext(args.logfile)[0]

    want = args.all or any([args.kml, args.gpx, args.geojson, args.csv, args.excel])
    if not want:
        ap.print_help()
        sys.exit(0)

    if args.all or args.kml:
        write_kml(rows,     args.kml     or stem + '.kml')
    if args.all or args.gpx:
        write_gpx(rows,     args.gpx     or stem + '.gpx')
    if args.all or args.geojson:
        write_geojson(rows, args.geojson or stem + '.geojson')
    if args.all or args.csv:
        write_csv(rows,     args.csv     or stem + '.csv')
    if args.all or args.excel:
        write_excel(rows,   args.excel   or stem + '.xlsx')

if __name__ == '__main__':
    main()
