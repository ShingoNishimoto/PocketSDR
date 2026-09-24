#!/usr/bin/env python3
"""drone_ppp_export.py - Export one or more offline drone PPP reprocessing
runs (rnx2rtkp -pppdbg <file> ... -o <pos>) into a single .xlsx workbook,
one sheet per run -- same per-epoch columns as gs_ppp_export.py's per-
ephemeris build() (position/quality/clock/PPP-KF diagnostics), but for the
drone's main (clock-fixed) solver only, no onboard reference-solver
comparison.

Usage:
    python3 drone_ppp_export.py work/20260710101607 \\
        --run Broadcast main_brdc_newclk.pos pppdbg_main_brdc_newclk.csv \\
        --run UltraRapid main_ultra_newclk.pos pppdbg_main_ultra_newclk.csv \\
        --run RapidPrecise main_precise_newclk.pos pppdbg_main_precise_newclk.csv \\
        -o main_newclk_ppp_results.xlsx
"""
import sys, os, argparse
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gs_ppp_export import build  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('session_dir')
    ap.add_argument('--run', nargs=3, action='append', metavar=('SHEET', 'POSFILE', 'PPPDBGFILE'),
                    required=True, help='sheet name + .pos + pppdbg .csv (relative to session_dir); repeatable')
    ap.add_argument('-o', '--out', default='main_newclk_ppp_results.xlsx')
    args = ap.parse_args()

    sd = args.session_dir.rstrip('/')
    out_path = os.path.join(sd, args.out)

    with pd.ExcelWriter(out_path, engine='openpyxl') as xw:
        for sheet, posfile, pppdbgfile in args.run:
            df = build(sd, sheet, posfile, pppdbgfile)
            df.to_excel(xw, sheet_name=sheet, index=False)
    print(f'wrote {out_path}')


if __name__ == '__main__':
    main()
