#!/usr/bin/env python3
#
#  pocket_rtk_compare.py - Diff onboard AOWR main/reference solver trajectories
#  against an independent RTKLIB kinematic-RTK ground truth (e.g. rover vs a
#  Canberra CORS station).
#
#  Double-difference RTK (rnx2rtkp against a CORS base) cancels the receiver
#  clock in the differencing, so a Canberra RTK run cannot itself show an
#  "AOWR-fixed-clock vs free-clock" difference -- there is no clock parameter
#  at that level. What it gives is a precise, independent ground-truth
#  trajectory. This script diffs the onboard
#    $POS    - main solver    (AOWR clock-fixed, clock shared from AOWR)
#    $REFPOS - reference solver (free clock, normal 4-unknown estimation)
#  against that RTK truth, epoch by epoch, to quantify how much (if any)
#  position accuracy is lost by trusting the AOWR-shared clock instead of
#  estimating it independently -- i.e. the actual time-transfer accuracy
#  check.
#
#  Workflow:
#    python3 pocket_obs2rnx.py pocket.log -o rover.obs
#    # download base RINEX for the flight window, e.g. GA CORS CANB:
#    #   https://gnss.ga.gov.au/rinex/
#    rnx2rtkp -k rtk_kinematic.conf rover.obs canb_base.obs canb_base.nav -o rtk.pos
#    python3 pocket_rtk_compare.py pocket.log rtk.pos --out compare.png --csv compare.csv
#
import sys, argparse, os
from datetime import datetime, timezone, timedelta
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pocket_export as pe
import pocket_pos_plot as pp

SOLVER_COLOR = {'main': 'tab:red', 'ref': 'tab:blue'}
SOLVER_LABEL = {'main': 'Main solver (AOWR clock-fixed)',
                'ref':  'Reference solver (free clock)'}


# ── time alignment ─────────────────────────────────────────────────────────────

def calendar_ts(r):
    """$POS/$REFPOS row's GPST calendar fields -> Unix-like timestamp.

    Labels the GPST calendar as UTC (no leap-second shift), matching
    pocket_pos_plot._gpst_to_ts's convention for RTKLIB .pos files, so the
    two absolute time axes line up directly.
    """
    si = int(r['sec'])
    frac = r['sec'] - si
    extra_min, si = divmod(si, 60)
    dt = datetime(r['year'], r['month'], r['day'], r['hour'], r['min'], 0,
                  tzinfo=timezone.utc) + timedelta(minutes=extra_min, seconds=si)
    return dt.timestamp() + frac


def interp_truth(truth_t, truth_xyz, t, max_gap):
    """Linearly interpolate truth ECEF (3,N) onto query times t (len M).

    Returns (xyz (3,M), valid (M,) bool) -- valid is False outside the truth
    time span or where the surrounding truth epochs are more than max_gap
    seconds apart (avoids silently interpolating across an RTK dropout).
    """
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


def rows_to_arrays(rows):
    t   = np.array([calendar_ts(r) for r in rows])
    xyz = np.array([[r['ecef_x'] for r in rows],
                     [r['ecef_y'] for r in rows],
                     [r['ecef_z'] for r in rows]])
    q   = np.array([r['q'] for r in rows], dtype=int)
    return t, xyz, q


def solver_vs_truth(t, xyz, truth_t, truth_xyz, lat0, lon0, max_gap):
    txyz, valid = interp_truth(truth_t, truth_xyz, t, max_gap)
    dx, dy, dz = xyz[0] - txyz[0], xyz[1] - txyz[1], xyz[2] - txyz[2]
    n, e, u = pp.ecef_to_neu(dx, dy, dz, lat0, lon0)
    return n, e, u, valid


def print_delta(main_stats, ref_stats):
    (mn, me, mu), (rn, re, ru) = main_stats, ref_stats
    m3d = np.sqrt(np.mean(mn**2 + me**2 + mu**2))
    r3d = np.sqrt(np.mean(rn**2 + re**2 + ru**2))
    print(f'\nMain solver (AOWR-fixed) 3D RMS vs RTK truth: {m3d:.3f} m')
    print(f'Reference solver (free clock) 3D RMS vs RTK truth: {r3d:.3f} m')
    print(f'Delta (main - ref): {m3d - r3d:+.3f} m '
          f'({"AOWR clock costs accuracy" if m3d > r3d else "AOWR clock is as good or better"})')


def main():
    ap = argparse.ArgumentParser(
        description='Diff onboard AOWR main/reference solver trajectories '
                    'against an RTKLIB kinematic-RTK ground truth')
    ap.add_argument('pocket_log', help='pocket.log (contains $POS, $REFPOS, $LOG SC_AOWR)')
    ap.add_argument('rtk_pos', help='RTKLIB rnx2rtkp .pos file (rover vs CORS base, e.g. CANB)')
    ap.add_argument('--truth-minq', type=int, default=2,
                    help='min RTK quality to trust as truth: 1=fix, 2=float (default 2)')
    ap.add_argument('--max-gap', type=float, default=2.0,
                    help='max seconds between truth epochs to interpolate across (default 2.0)')
    ap.add_argument('--csv', default=None, help='write per-epoch comparison CSV')
    ap.add_argument('--out', default=None, help='save figure to PNG instead of displaying')
    args = ap.parse_args()

    main_rows = pe.parse_log(args.pocket_log)
    ref_rows  = pe.parse_refpos_log(args.pocket_log)
    aowr_hist = pe.parse_sc_aowr_log(args.pocket_log)

    if not main_rows:
        print('No $POS (main solver) records found.', file=sys.stderr)
        sys.exit(1)
    if not ref_rows:
        print('Warning: no $REFPOS (reference solver) records found -- was '
              '-ps_sc active? Only the main-solver-vs-truth comparison will '
              'be produced.', file=sys.stderr)

    truth = pp._parse_rtklib_pos(args.rtk_pos)
    truth = truth[(truth['q'] >= 1) & (truth['q'] <= args.truth_minq)]
    if len(truth) == 0:
        print(f'No truth epochs in {args.rtk_pos} meet --truth-minq {args.truth_minq}.',
              file=sys.stderr)
        sys.exit(1)

    lat0 = float(np.mean(truth['lat']))
    lon0 = float(np.mean(truth['lon']))
    hgt0 = float(np.mean(truth['hgt']))
    truth_xyz = np.array(pp.llh_to_ecef(truth['lat'], truth['lon'], truth['hgt']))

    main_t, main_xyz, main_q = rows_to_arrays(main_rows)
    mn, me, mu, mvalid = solver_vs_truth(main_t, main_xyz, truth['t'], truth_xyz,
                                          lat0, lon0, args.max_gap)
    print(f'Main solver: {len(main_rows)} epochs, {int(np.sum(mvalid))} '
          f'matched to truth (within {args.max_gap}s)')
    if np.any(mvalid):
        pp.print_stats(SOLVER_LABEL['main'], mn[mvalid], me[mvalid], mu[mvalid])

    ref_t = ref_xyz = ref_q = None
    rn = re = ru = rvalid = None
    if ref_rows:
        ref_t, ref_xyz, ref_q = rows_to_arrays(ref_rows)
        rn, re, ru, rvalid = solver_vs_truth(ref_t, ref_xyz, truth['t'], truth_xyz,
                                              lat0, lon0, args.max_gap)
        print(f'\nReference solver: {len(ref_rows)} epochs, {int(np.sum(rvalid))} '
              f'matched to truth (within {args.max_gap}s)')
        if np.any(rvalid):
            pp.print_stats(SOLVER_LABEL['ref'], rn[rvalid], re[rvalid], ru[rvalid])

        if np.any(mvalid) and np.any(rvalid):
            print_delta((mn[mvalid], me[mvalid], mu[mvalid]),
                        (rn[rvalid], re[rvalid], ru[rvalid]))

    # ── CSV ──────────────────────────────────────────────────────────────────
    if args.csv:
        import csv as csvmod
        with open(args.csv, 'w', newline='') as f:
            w = csvmod.writer(f)
            w.writerow(['t_abs', 'datetime_gpst', 'solver', 'q', 'valid',
                        'n_err_m', 'e_err_m', 'u_err_m', 'sc_aowr_clk_s',
                        'sc_aowr_drift'])
            for i, r in enumerate(main_rows):
                aowr = aowr_hist.get(r['t'], {})
                w.writerow([f'{main_t[i]:.3f}', pe.gpst_str(r), 'main', main_q[i],
                            int(mvalid[i]), f'{mn[i]:.4f}', f'{me[i]:.4f}',
                            f'{mu[i]:.4f}', aowr.get('clk', ''), aowr.get('drift', '')])
            if ref_rows:
                for i, r in enumerate(ref_rows):
                    w.writerow([f'{ref_t[i]:.3f}', pe.gpst_str(r), 'ref', ref_q[i],
                                int(rvalid[i]), f'{rn[i]:.4f}', f'{re[i]:.4f}',
                                f'{ru[i]:.4f}', '', ''])
        print(f'\nCSV: {args.csv}')

    # ── plot ─────────────────────────────────────────────────────────────────
    t0 = main_t[0]
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle('AOWR time-transfer check: main vs reference solver vs RTK truth',
                 fontsize=11)
    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.5, wspace=0.3,
                           height_ratios=[1, 1, 1, 1])
    ax_n   = fig.add_subplot(gs[0, 0])
    ax_e   = fig.add_subplot(gs[1, 0])
    ax_u   = fig.add_subplot(gs[2, 0])
    ax_clk = fig.add_subplot(gs[3, 0])
    ax_ne  = fig.add_subplot(gs[0:2, 1], aspect='equal')
    ax_3d  = fig.add_subplot(gs[2:4, 1])

    def plot_series(ax, t, val, valid, tag, ylabel):
        tm = (t - t0) / 60.0
        ax.scatter(tm[valid], val[valid], s=4, c=SOLVER_COLOR[tag],
                   label=SOLVER_LABEL[tag])
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)

    plot_series(ax_n, main_t, mn, mvalid, 'main', 'N err (m)')
    plot_series(ax_e, main_t, me, mvalid, 'main', 'E err (m)')
    plot_series(ax_u, main_t, mu, mvalid, 'main', 'U err (m)')
    ax_ne.scatter(me[mvalid], mn[mvalid], s=5, c=SOLVER_COLOR['main'],
                 label=SOLVER_LABEL['main'])
    if ref_rows:
        plot_series(ax_n, ref_t, rn, rvalid, 'ref', 'N err (m)')
        plot_series(ax_e, ref_t, re, rvalid, 'ref', 'E err (m)')
        plot_series(ax_u, ref_t, ru, rvalid, 'ref', 'U err (m)')
        ax_ne.scatter(re[rvalid], rn[rvalid], s=5, c=SOLVER_COLOR['ref'],
                     label=SOLVER_LABEL['ref'])
    ax_n.legend(fontsize=7); ax_u.set_xlabel('minutes from session start')
    ax_ne.set_xlabel('E err (m)'); ax_ne.set_ylabel('N err (m)')
    ax_ne.grid(alpha=0.3); ax_ne.legend(fontsize=7)

    # AOWR clock offset history (as actually applied to the main solver)
    aowr_t = np.array([r['t'] for r in main_rows])
    aowr_clk = np.array([aowr_hist.get(r['t'], {}).get('clk', np.nan) for r in main_rows])
    ax_clk.plot((aowr_t - aowr_t[0]) / 60.0, aowr_clk * 1e9, '.', ms=3, color='tab:green')
    ax_clk.set_ylabel('AOWR clk (ns)'); ax_clk.set_xlabel('minutes from session start')
    ax_clk.grid(alpha=0.3)

    # Running 3D error RMS (cumulative) for a quick "is this converging" check
    if np.any(mvalid):
        m3d_run = np.sqrt(np.cumsum((mn[mvalid]**2 + me[mvalid]**2 + mu[mvalid]**2)) /
                          (np.arange(int(np.sum(mvalid))) + 1))
        ax_3d.plot((main_t[mvalid] - t0) / 60.0, m3d_run, color=SOLVER_COLOR['main'],
                  label=SOLVER_LABEL['main'])
    if ref_rows and np.any(rvalid):
        r3d_run = np.sqrt(np.cumsum((rn[rvalid]**2 + re[rvalid]**2 + ru[rvalid]**2)) /
                          (np.arange(int(np.sum(rvalid))) + 1))
        ax_3d.plot((ref_t[rvalid] - t0) / 60.0, r3d_run, color=SOLVER_COLOR['ref'],
                  label=SOLVER_LABEL['ref'])
    ax_3d.set_ylabel('cumulative 3D RMS (m)'); ax_3d.set_xlabel('minutes from session start')
    ax_3d.grid(alpha=0.3); ax_3d.legend(fontsize=7)

    if args.out:
        fig.savefig(args.out, dpi=150)
        print(f'\nFigure: {args.out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
