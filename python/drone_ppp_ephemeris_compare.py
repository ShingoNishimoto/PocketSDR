#!/usr/bin/env python3
"""drone_ppp_ephemeris_compare.py - Compare drone-side clock-fixed PPP
kinematic solutions (rnx2rtkp -clkfix, one run per ephemeris type) against
an independent RTK ground truth (rover vs CORS base), epoch by epoch.

Pass a pppdbg file with each --run to get the same rich per-epoch columns
as drone_ppp_export.py (position/quality/clock/PPP-KF diagnostics), with
n_err_m/e_err_m/u_err_m/valid_vs_truth merged in on top -- this is what
--excel/--csv write out. Without a pppdbg file, --run falls back to a
minimal position/quality-only frame (still gets the error columns, just
none of the PPP-KF diagnostics).

Usage:
    python3 drone_ppp_ephemeris_compare.py work/20260710101607 rtk.pos \\
        --run broadcast main_brdc_newclk.pos pppdbg_main_brdc_newclk.csv \\
        --run ultra-rapid main_ultra_newclk.pos pppdbg_main_ultra_newclk.csv \\
        --run rapid-precise main_precise_newclk.pos pppdbg_main_precise_newclk.csv \\
        --out eph_compare.png --excel eph_compare.xlsx
"""
import sys, os, argparse
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pocket_pos_plot as pp
from gs_ppp_export import build as build_rich  # noqa: E402


def colors_for(n):
    """n distinct colors -- tab10 cycles/repeats past 10 runs, which is a
    saner failure mode than zip()-truncating and silently dropping runs."""
    cmap = plt.get_cmap('tab10')
    return [cmap(i % 10) for i in range(n)]


def interp_truth(truth_t, truth_xyz, t, max_gap):
    order = np.argsort(truth_t)
    tt, xyz = truth_t[order], truth_xyz[:, order]
    out = np.full((3, len(t)), np.nan)
    if len(tt) < 2:
        return out, np.zeros(len(t), dtype=bool)
    valid = (t >= tt[0]) & (t <= tt[-1])
    idx = np.clip(np.searchsorted(tt, t), 1, len(tt) - 1)
    valid &= (tt[idx] - tt[idx - 1]) <= max_gap
    for k in range(3):
        out[k, valid] = np.interp(t[valid], tt, xyz[k])
    return out, valid


def datetime_gpst_to_ts(s):
    """Same convention as pocket_pos_plot._gpst_to_ts: GPST calendar labeled
    as UTC, no leap-second shift -- must match so t_abs lines up with the
    truth .pos file's own timestamps."""
    return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f').replace(tzinfo=timezone.utc).timestamp()


def load_clkcmp(path):
    """Parse a -clkcmp companion CSV (# time,aowr_dtr_s,ppp_dtr_s,diff_s,
    diff_m) -- aowr_dtr_s is the actual shared/external clock value that
    was fed into -clkfix for that epoch (from lookup_ext_clk() in
    postpos.c), as distinct from the PPP-estimated dtr_s/ppp_clk_m already
    in the rich columns. Only written for epochs where the external clock
    lookup succeeded (e.g. absent during a staged pre-clock-fix window --
    that's expected, not a bug). Combined-mode (forward+backward) runs
    write each epoch twice; dedup to the first (forward-pass) occurrence,
    same convention as gs_ppp_export.load_dedup_csv for pppdbg.
    Returns {t_abs: {'aowr_dtr_s':, 'ppp_dtr_s_clkcmp':, 'clkcmp_diff_s':,
    'clkcmp_diff_m':}}."""
    rows = {}
    with open(path) as f:
        next(f)  # header comment
        for line in f:
            p = line.strip().split(',')
            if len(p) != 5:
                continue
            if p[0] in rows:
                continue
            d, t = p[0].split(' ')
            rows[p[0]] = dict(
                t_abs=pp._gpst_to_ts(d, t), aowr_dtr_s=float(p[1]),
                ppp_dtr_s_clkcmp=float(p[2]), clkcmp_diff_s=float(p[3]),
                clkcmp_diff_m=float(p[4]))
    return rows


def load_run(sd, label, posfile, pppdbgfile, clkcmpfile=None):
    """Returns a DataFrame with a 't_abs'/'lat_deg'/'lon_deg'/'hgt_m'/
    'quality' column at minimum -- the full drone_ppp_export.py column set
    when pppdbgfile is given, a minimal position/quality frame otherwise.
    clkcmpfile (only meaningful with pppdbgfile) adds aowr_dtr_s and the
    PPP-vs-AOWR diff columns from load_clkcmp(), matched by t_abs."""
    if pppdbgfile:
        df = build_rich(sd, label, posfile, pppdbgfile)
        df['t_abs'] = df['datetime_gpst'].apply(datetime_gpst_to_ts)
        if clkcmpfile:
            clkcmp = load_clkcmp(os.path.join(sd, clkcmpfile))
            by_t = {round(v['t_abs']): v for v in clkcmp.values()}
            cols = ['aowr_dtr_s', 'ppp_dtr_s_clkcmp', 'clkcmp_diff_s', 'clkcmp_diff_m']
            for c in cols:
                df[c] = [by_t.get(round(t), {}).get(c, float('nan')) for t in df['t_abs']]
            n_matched = df['aowr_dtr_s'].notna().sum()
            print(f'{label}: matched {n_matched}/{len(df)} rows to {clkcmpfile}')
    else:
        data = pp._parse_rtklib_pos(os.path.join(sd, posfile))
        df = pd.DataFrame({
            't_abs': data['t'], 'quality': data['q'], 'n_sats': data['ns'],
            'lat_deg': data['lat'], 'lon_deg': data['lon'], 'hgt_m': data['hgt'],
        })
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('session_dir')
    ap.add_argument('truth_pos', help='RTK ground-truth .pos filename, relative to session_dir')
    ap.add_argument('--run', nargs='+', action='append', metavar='LABEL POSFILE [PPPDBGFILE [CLKCMPFILE]]',
                    required=True, help='label + .pos filename + optional pppdbg csv + '
                         'optional -clkcmp csv (all relative to session_dir); repeatable. '
                         'Include the pppdbg file to get the full rich column set in '
                         '--excel/--csv; include the clkcmp file too (needs pppdbg) to add '
                         'aowr_dtr_s (the actual shared/external clock fed into -clkfix) '
                         'alongside the PPP-estimated dtr_s/ppp_clk_m columns.')
    ap.add_argument('--truth-minq', type=int, default=2)
    ap.add_argument('--max-gap', type=float, default=2.0)
    ap.add_argument('--csv', default=None)
    ap.add_argument('--excel', default=None,
                    help='write .xlsx with one sheet per --run label plus a '
                         'Summary sheet, instead of one long CSV -- easier to '
                         'compare runs side by side')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    sd = args.session_dir.rstrip('/')
    truth = pp._parse_rtklib_pos(os.path.join(sd, args.truth_pos))
    truth = truth[(truth['q'] >= 1) & (truth['q'] <= args.truth_minq)]
    if len(truth) == 0:
        sys.exit(f'no truth epochs in {args.truth_pos} meet --truth-minq {args.truth_minq}')

    lat0, lon0, hgt0 = float(np.mean(truth['lat'])), float(np.mean(truth['lon'])), float(np.mean(truth['hgt']))
    truth_xyz = np.array(pp.llh_to_ecef(truth['lat'], truth['lon'], truth['hgt']))

    results = {}
    for run_args in args.run:
        if len(run_args) == 2:
            label, posfile, pppdbgfile, clkcmpfile = *run_args, None, None
        elif len(run_args) == 3:
            label, posfile, pppdbgfile = run_args
            clkcmpfile = None
        elif len(run_args) == 4:
            label, posfile, pppdbgfile, clkcmpfile = run_args
        else:
            sys.exit(f'--run expects LABEL POSFILE [PPPDBGFILE [CLKCMPFILE]], got {run_args}')

        df = load_run(sd, label, posfile, pppdbgfile, clkcmpfile)
        if len(df) == 0:
            print(f'{label}: no epochs in {posfile}, skipping', file=sys.stderr)
            continue
        t_abs = df['t_abs'].to_numpy()
        xyz = np.array(pp.llh_to_ecef(df['lat_deg'].to_numpy(), df['lon_deg'].to_numpy(), df['hgt_m'].to_numpy()))
        txyz, valid = interp_truth(truth['t'], truth_xyz, t_abs, args.max_gap)
        dx, dy, dz = xyz[0] - txyz[0], xyz[1] - txyz[1], xyz[2] - txyz[2]
        n, e, u = pp.ecef_to_neu(dx, dy, dz, lat0, lon0)
        df['n_err_m'], df['e_err_m'], df['u_err_m'], df['valid_vs_truth'] = n, e, u, valid

        print(f'{label}: {len(df)} epochs, {int(np.sum(valid))} matched to truth '
              f'(within {args.max_gap}s)')
        if np.any(valid):
            pp.print_stats(label, n[valid], e[valid], u[valid])
        q = df['quality'].to_numpy() if 'quality' in df else np.full(len(df), np.nan)
        results[label] = dict(df=df, t=t_abs, n=n, e=e, u=u, valid=valid, q=q)

    if not results:
        sys.exit('no runs produced any matched epochs')

    # ── summary table (RMS 3D + mean bias per axis) ─────────────────────────
    summary_rows = []
    print('\n=== Summary: 3D RMS vs RTK truth ===')
    for label, r in results.items():
        v = r['valid']
        if not np.any(v):
            continue
        n, e, u = r['n'][v], r['e'][v], r['u'][v]
        rms3d = np.sqrt(np.mean(n**2 + e**2 + u**2))
        print(f'  {label:16s} n={int(np.sum(v)):5d}  3D RMS = {rms3d:7.3f} m')
        summary_rows.append(dict(
            label=label, n_epochs=int(np.sum(v)),
            n_bias_m=np.mean(n), e_bias_m=np.mean(e), u_bias_m=np.mean(u),
            n_rms_m=np.sqrt(np.mean(n**2)), e_rms_m=np.sqrt(np.mean(e**2)),
            u_rms_m=np.sqrt(np.mean(u**2)), rms_3d_m=rms3d))

    # ── CSV (long format: all runs concatenated, distinguished by 'label') ──
    if args.csv:
        frames = []
        for label, r in results.items():
            df = r['df'].copy()
            df.insert(0, 'label', label)
            frames.append(df)
        pd.concat(frames, ignore_index=True).to_csv(os.path.join(sd, args.csv), index=False)
        print(f'\nCSV: {args.csv}')

    # ── Excel (one sheet per run, full column set, + a Summary sheet) ──────
    if args.excel:
        excel_path = os.path.join(sd, args.excel)
        with pd.ExcelWriter(excel_path, engine='openpyxl') as xw:
            pd.DataFrame(summary_rows).to_excel(xw, sheet_name='Summary', index=False)
            for label, r in results.items():
                sheet = label[:31]  # Excel sheet-name length limit
                r['df'].to_excel(xw, sheet_name=sheet, index=False)
        print(f'Excel: {excel_path}')

    # ── plot ─────────────────────────────────────────────────────────────────
    t0 = min(r['t'][0] for r in results.values())
    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    ax_n, ax_e, ax_u, ax_3d = axes
    for (label, r), c in zip(results.items(), colors_for(len(results))):
        v = r['valid']
        tm = (r['t'] - t0) / 60.0
        ax_n.scatter(tm[v], r['n'][v], s=5, c=c, label=label)
        ax_e.scatter(tm[v], r['e'][v], s=5, c=c, label=label)
        ax_u.scatter(tm[v], r['u'][v], s=5, c=c, label=label)
        if np.any(v):
            run3d = np.sqrt(np.cumsum(r['n'][v]**2 + r['e'][v]**2 + r['u'][v]**2) /
                            (np.arange(int(np.sum(v))) + 1))
            ax_3d.plot(tm[v], run3d, c=c, label=label)
    ax_n.set_ylabel('N err (m)'); ax_n.grid(alpha=0.3); ax_n.legend(fontsize=8)
    ax_e.set_ylabel('E err (m)'); ax_e.grid(alpha=0.3)
    ax_u.set_ylabel('U err (m)'); ax_u.grid(alpha=0.3)
    ax_3d.set_ylabel('cumulative 3D RMS (m)'); ax_3d.set_xlabel('minutes from session start')
    ax_3d.grid(alpha=0.3)
    fig.suptitle(f'{sd}: drone clock-fixed PPP vs RTK truth, by ephemeris type')
    fig.tight_layout()

    if args.out:
        out_path = os.path.join(sd, args.out)
        fig.savefig(out_path, dpi=150)
        print(f'\nFigure: {out_path}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
