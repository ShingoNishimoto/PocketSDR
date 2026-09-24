#!/usr/bin/env python3
"""Export a GS (ground-station) offline PPP reprocessing run (rapid precise
or IGS ultra-rapid ephemeris, produced via rnx2rtkp -pppdbg) into a
pocket.xlsx-equivalent workbook, alongside the session's own real-time
broadcast-ephemeris run, for direct comparison.

Works for both GS processing modes:
  - ppp-fixed (position pinned to a known survey point, only clock/trop/amb
    estimated): pass --prefix gs
  - ppp-static (position estimated via random walk from its own SPP seed):
    pass --prefix gs_static

Expects these files under SESSION_DIR (produced by rnx2rtkp -pppdbg
<file> -k <conf> rover.obs rover.nav <sp3> [<clk>] -o <prefix>_precise.pos,
run twice with rapid-precise and ultra-rapid ephemeris):
    <prefix>_precise.pos, pppdbg_<prefix>_precise.csv
    <prefix>_ultra.pos,   pppdbg_<prefix>_ultra.csv
    <brdc_xlsx> (the session's own real-time export, e.g. pocket_gs.xlsx or
    pocket.xlsx -- filename varies by session, sheet name is normally 'PPP')

Writes:
    <prefix>_precise_ppp.csv, <prefix>_ultra_ppp.csv
    <prefix>_ppp_results.xlsx  (Broadcast/RapidPrecise/UltraRapid sheets)
    <prefix>_clock_compare.png (clock offset + PPP clock std, common window)

Usage:
    python3 gs_ppp_export.py <session_dir> <prefix> <brdc_xlsx> [brdc_sheet]

Example (GS fixed-position mode, session work/gs/3-5):
    python3 gs_ppp_export.py ../work/gs/3-5 gs pocket_gs.xlsx

Example (GS static/position-estimated mode, session work/gs/3-4):
    python3 gs_ppp_export.py ../work/gs/3-4 gs_static pocket.xlsx

Background: clk_m in pppdbg is the PPP KF's POSTERIOR clock state
(x[IC(0)] in RTKLIB_fork's ppp.c), not necessarily the raw input clock --
for -clkfix runs the two can differ under the soft-constraint model. This
export always reports the posterior (dtr_s / ppp_clk_m), consistent with
pocket_export.py's own column semantics.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from pocket_export import parse_rtklib_pos_as_ref, make_df  # noqa: E402

CLIGHT = 299792458.0

PPPDBG_TO_ROW = {
    'kf_pos_std_m': 'ppp_pos_std', 'ndual': 'ppp_ndual', 'clk_std_m': 'ppp_clk_std',
    'zwd_m': 'ppp_zwd', 'zwd_std_m': 'ppp_zwd_std', 'amb_std_m': 'ppp_amb_std',
    'nconv': 'ppp_nconv', 'res_c_m': 'ppp_res_c',
    'gdop': 'gdop', 'pdop': 'pdop', 'hdop': 'hdop', 'vdop': 'vdop',
}


def load_dedup_csv(path, n_expect):
    """Parse a -pppdbg/-clkcmp style CSV, keeping only the first (forward-
    pass) occurrence of each timestamp -- combres() in RTKLIB_fork's
    postpos.c carries the forward pass's clock/KF state through to the
    combined solution whenever forward and backward cover the same epoch."""
    rows = {}
    with open(path) as f:
        header = next(f).lstrip('# \n').strip().split(',')
        for line in f:
            p = line.strip().split(',')
            if len(p) != n_expect:
                continue
            t = p[0]
            if t in rows:
                continue
            rows[t] = dict(zip(header[1:], (float(x) for x in p[1:])))
    return rows


def row_time_key(r):
    return (f"{r['year']:04d}/{r['month']:02d}/{r['day']:02d} "
            f"{r['hour']:02d}:{r['min']:02d}:{r['sec']:06.3f}")


def merge_pppdbg(rows, pppdbg):
    n_matched = 0
    for r in rows:
        d = pppdbg.get(row_time_key(r))
        if d is None:
            continue
        r['dtr'] = d['clk_m'] / CLIGHT
        r['ppp_clk'] = d['clk_m']
        for src, dst in PPPDBG_TO_ROW.items():
            r[dst] = d[src]
        n_matched += 1
    return n_matched


def add_drift(rows):
    """Finite-difference dtr_drift_s -- PPP's clock is a white-noise KF
    state re-seeded each epoch, not a modeled drift-rate state, so this is
    a derived post-hoc quantity. NaN at the first row / gaps >5s."""
    rows.sort(key=lambda r: r['t'])
    for i, r in enumerate(rows):
        r['dtr_drift'] = float('nan')
        if i == 0:
            continue
        prev = rows[i - 1]
        dt = r['t'] - prev['t']
        if 0 < dt <= 5.0 and 'dtr' in r and 'dtr' in prev:
            r['dtr_drift'] = (r['dtr'] - prev['dtr']) / dt


def build(session_dir, label, pos_file, pppdbg_file):
    rows = parse_rtklib_pos_as_ref(f'{session_dir}/{pos_file}')
    pppdbg = load_dedup_csv(f'{session_dir}/{pppdbg_file}', 15)
    n = merge_pppdbg(rows, pppdbg)
    print(f'{label}: matched {n}/{len(rows)} rows to {pppdbg_file}')
    add_drift(rows)
    return make_df(rows)


def hms_to_s(t):
    h, m, s = t.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)


def clean(df, dtr_col, std_col):
    d = df[['datetime_gpst', dtr_col, std_col]].dropna().copy()
    tstr = d['datetime_gpst'].astype(str).str.split('T').str[1]
    d['t_s'] = tstr.apply(hms_to_s)
    return d.sort_values('t_s')


def tail_stats(label, d, dtr_col='dtr_s', std_col='ppp_clk_std_m'):
    n = len(d)
    tail = d.iloc[n // 3:]
    print(f'{label:32s} n={n:5d}  dtr std(tail)={tail[dtr_col].std()*1e3:8.4f} ms  '
          f'mean clk_std(tail)={tail[std_col].mean():7.4f} m')


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    session_dir, prefix, brdc_xlsx = sys.argv[1], sys.argv[2], sys.argv[3]
    brdc_sheet = sys.argv[4] if len(sys.argv) > 4 else 'PPP'

    df_precise = build(session_dir, f'{prefix}-precise', f'{prefix}_precise.pos',
                        f'pppdbg_{prefix}_precise.csv')
    df_ultra = build(session_dir, f'{prefix}-ultra', f'{prefix}_ultra.pos',
                      f'pppdbg_{prefix}_ultra.csv')

    df_precise.to_csv(f'{session_dir}/{prefix}_precise_ppp.csv', index=False, float_format='%.9g')
    df_ultra.to_csv(f'{session_dir}/{prefix}_ultra_ppp.csv', index=False, float_format='%.9g')

    brdc = pd.read_excel(f'{session_dir}/{brdc_xlsx}', sheet_name=brdc_sheet)

    xlsx_path = f'{session_dir}/{prefix}_ppp_results.xlsx'
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as xw:
        brdc.to_excel(xw, sheet_name='Broadcast', index=False)
        df_precise.to_excel(xw, sheet_name='RapidPrecise', index=False)
        df_ultra.to_excel(xw, sheet_name='UltraRapid', index=False)
    print(f'wrote {xlsx_path}')

    b = clean(brdc, 'dtr_s', 'ppp_clk_std_m')
    p = clean(df_precise, 'dtr_s', 'ppp_clk_std_m')
    u = clean(df_ultra, 'dtr_s', 'ppp_clk_std_m')
    t0 = min(b['t_s'].min(), p['t_s'].min(), u['t_s'].min())

    # Clip all series to the shortest common window: the precise/ultra-rapid
    # offline reprocessing runs can stop earlier than the broadcast baseline
    # (observed once as a numerical-convergence limit in RTKLIB_fork's PPP
    # solver, unrelated to ephemeris choice -- confirmed by reproducing the
    # same cutoff with broadcast ephemeris under the same strict config).
    # Without this, the longer series would visually/statistically run past
    # where the others stop, which isn't a fair same-period comparison.
    t_end = min(b['t_s'].max(), p['t_s'].max(), u['t_s'].max())
    print(f'clipping comparison plot to common window: t0={t0:.1f} t_end={t_end:.1f} '
          f'({t_end - t0:.1f}s, vs broadcast\'s own {b["t_s"].max() - t0:.1f}s)')
    b = b[b['t_s'] <= t_end]
    p = p[p['t_s'] <= t_end]
    u = u[u['t_s'] <= t_end]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    ax = axes[0]
    ax.plot(b['t_s'] - t0, b['dtr_s'] * 1e3, '.', ms=3, label='broadcast eph (real-time)')
    ax.plot(p['t_s'] - t0, p['dtr_s'] * 1e3, '.', ms=3, label='rapid precise (GFZ MGX)')
    ax.plot(u['t_s'] - t0, u['dtr_s'] * 1e3, '.', ms=3, label='ultra-rapid (IGU predicted + QZS brdc)')
    ax.set_ylabel('GS receiver clock offset (ms)')
    ax.set_title(f'{session_dir} [{prefix}]: broadcast vs rapid vs ultra-rapid ephemeris (common window)')
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(b['t_s'] - t0, b['ppp_clk_std_m'], '.', ms=3, label='broadcast eph')
    ax.plot(p['t_s'] - t0, p['ppp_clk_std_m'], '.', ms=3, label='rapid precise')
    ax.plot(u['t_s'] - t0, u['ppp_clk_std_m'], '.', ms=3, label='ultra-rapid')
    ax.set_ylabel('PPP clock std (m)')
    ax.set_xlabel('Elapsed time (s)')
    ax.set_ylim(0, 2)
    ax.legend(fontsize=8)

    plt.tight_layout()
    out = f'{session_dir}/{prefix}_clock_compare.png'
    plt.savefig(out, dpi=120)
    print(f'wrote {out}')

    print()
    tail_stats('broadcast (real-time, clipped)', b)
    tail_stats('rapid precise (GFZ MGX)', p)
    tail_stats('ultra-rapid (IGU+QZS)', u)


if __name__ == '__main__':
    main()
