#!/usr/bin/env python3
#
#  pocket_rtk_export.py - Export an RTKLIB rnx2rtkp .pos file's full
#  position/state history to CSV and/or Excel.
#
#  Complements pocket_rtk_compare.py's --csv (which only holds the diff
#  against the onboard main/reference solvers): this exports the raw RTK
#  solution itself -- position, quality, satellite count, formal sigmas
#  (including the covariance cross-terms RTKLIB reports, sdne/sdeu/sdun),
#  differential age and AR ratio -- for every epoch.
#
#  Requires the .pos file to have been produced with out-timeform=hms (see
#  conf/rtk_kinematic.conf and conf/rtk_static.conf) -- out-timeform=tow
#  (week/tow) columns are not parsed here, same requirement as
#  pocket_pos_plot.py / pocket_rtk_compare.py.
#
#  Usage:
#    python3 pocket_rtk_export.py rtk.pos --csv rtk.csv --excel rtk.xlsx
#    python3 pocket_rtk_export.py rtk.pos --all   (names both from rtk.pos stem)
#
import sys, argparse, os
from datetime import datetime, timezone
import numpy as np
import pandas as pd

SOLQ_LABEL = {0: 'NONE', 1: 'FIX', 2: 'FLOAT', 3: 'SBAS', 4: 'DGPS',
              5: 'SINGLE', 6: 'PPP', 7: 'DR'}

_WGS84_A  = 6378137.0
_WGS84_E2 = 2 / 298.257223563 - (1 / 298.257223563) ** 2

_FIELDS = ['lat', 'lon', 'hgt', 'q', 'ns', 'sdn', 'sde', 'sdu',
           'sdne', 'sdeu', 'sdun', 'age', 'ratio']
_INT_FIELDS = {'q', 'ns'}


GPS_EPOCH_TS = datetime(1980, 1, 6, tzinfo=timezone.utc).timestamp()


def gpst_week_tow(date_str, time_str):
    """rtk.pos 'YYYY/MM/DD' + 'HH:MM:SS.sss' (out-timesys=gpst) -> (week, tow_s)."""
    dt = datetime.strptime(f'{date_str} {time_str}', '%Y/%m/%d %H:%M:%S.%f')
    total = dt.replace(tzinfo=timezone.utc).timestamp() - GPS_EPOCH_TS
    week = int(total // 604800)
    return week, total - week * 604800


def llh_to_ecef(lat_deg, lon_deg, hgt_m):
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    N = _WGS84_A / np.sqrt(1.0 - _WGS84_E2 * np.sin(lat) ** 2)
    X = (N + hgt_m) * np.cos(lat) * np.cos(lon)
    Y = (N + hgt_m) * np.cos(lat) * np.sin(lon)
    Z = (N * (1.0 - _WGS84_E2) + hgt_m) * np.sin(lat)
    return X, Y, Z


def parse_rtk_pos(path):
    """Parse every column of an RTKLIB LLH .pos file (out-timeform=hms).

    Columns: date time lat lon hgt Q ns sdn sde sdu sdne sdeu sdun age ratio
    """
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%') or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 15:
                continue
            try:
                vals = [float(x) for x in p[2:15]]
            except ValueError:
                continue
            r = dict(zip(_FIELDS, vals))
            for k in _INT_FIELDS:
                r[k] = int(r[k])
            r['date'], r['time'] = p[0], p[1]
            rows.append(r)
    return rows


def make_df(rows):
    lat = np.array([r['lat'] for r in rows])
    lon = np.array([r['lon'] for r in rows])
    hgt = np.array([r['hgt'] for r in rows])
    X, Y, Z = llh_to_ecef(lat, lon, hgt)
    return pd.DataFrame([{
        'datetime_gpst': f"{r['date']}T{r['time']}",
        'gps_week': (_wt := gpst_week_tow(r['date'], r['time']))[0],
        'tow_s': _wt[1],
        'lat_deg': r['lat'], 'lon_deg': r['lon'], 'hgt_m': r['hgt'],
        'ecef_x_m': x, 'ecef_y_m': y, 'ecef_z_m': z,
        'quality': r['q'], 'quality_label': SOLQ_LABEL.get(r['q'], str(r['q'])),
        'n_sats': r['ns'],
        'sdn_m': r['sdn'], 'sde_m': r['sde'], 'sdu_m': r['sdu'],
        'sdne_m': r['sdne'], 'sdeu_m': r['sdeu'], 'sdun_m': r['sdun'],
        'age_s': r['age'], 'ar_ratio': r['ratio'],
    } for r, x, y, z in zip(rows, X, Y, Z)])


def write_csv(df, path):
    df.to_csv(path, index=False, float_format='%.9g')
    print(f'CSV:   {path}  ({len(df)} rows)')


def write_excel(df, path):
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        df.to_excel(xw, sheet_name='All', index=False)
        for label, q in [('FIX', 1), ('FLOAT', 2), ('DGPS', 4), ('SINGLE', 5)]:
            sub = df[df['quality'] == q]
            if len(sub):
                sub.to_excel(xw, sheet_name=label, index=False)
    print(f'Excel: {path}  ({len(df)} rows)')


def main():
    ap = argparse.ArgumentParser(
        description="Export an RTKLIB .pos file's position/state history to "
                    "CSV/Excel")
    ap.add_argument('posfile', help='RTKLIB rnx2rtkp .pos file (out-timeform=hms)')
    ap.add_argument('--csv', default=None, help='output CSV path')
    ap.add_argument('--excel', default=None, help='output Excel (.xlsx) path')
    ap.add_argument('--all', action='store_true',
                    help='write both CSV and Excel, named from posfile stem')
    args = ap.parse_args()

    rows = parse_rtk_pos(args.posfile)
    if not rows:
        print(f'No usable rows in {args.posfile} -- this needs out-timeform=hms '
              '(calendar date/time columns), not the week/tow format.',
              file=sys.stderr)
        sys.exit(1)
    df = make_df(rows)

    stem = os.path.splitext(args.posfile)[0]
    csv_path = args.csv or (f'{stem}.csv' if args.all else None)
    excel_path = args.excel or (f'{stem}.xlsx' if args.all else None)
    if not csv_path and not excel_path:
        excel_path = f'{stem}.xlsx'

    if csv_path:
        write_csv(df, csv_path)
    if excel_path:
        write_excel(df, excel_path)

    print(f'\n{len(df)} epochs:')
    for q, n in df['quality'].value_counts().sort_index().items():
        print(f'  {SOLQ_LABEL.get(q, str(q))}: {n} ({100 * n / len(df):.1f}%)')


if __name__ == '__main__':
    main()
