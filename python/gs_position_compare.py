#!/usr/bin/env python3
"""Compare a GS static-PPP run (position estimated, not pinned) against the
known surveyed position, for both rapid-precise and ultra-rapid ephemeris.
Only meaningful for ppp-static mode -- ppp-fixed pins position exactly (by
construction, 0 error every epoch), so there's nothing to score there.

Usage:
    python3 gs_position_compare.py <session_dir> <lat_deg> <lon_deg> <hgt_m> [prefix]

prefix defaults to gs_static, matching gs_ppp_export.py's --prefix gs_static
output naming (<prefix>_precise.pos / <prefix>_ultra.pos).

Example (truth LLH taken from the session's own -fixpos argument used for
the real-time ppp-fixed run, e.g. in real-time_ppp_l1l2_hybrid_station_fixpos.sh):
    python3 gs_position_compare.py ../work/gs/3-5 -35.281349827 149.116038399 581.281
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_pos(fn):
    rows = []
    with open(fn) as f:
        for line in f:
            if line.startswith('%') or not line.strip():
                continue
            p = line.split()
            t = p[0] + ' ' + p[1]
            rows.append((t, float(p[2]), float(p[3]), float(p[4])))
    return rows


def diffs(rows, truth_llh):
    la0, lo0, h0 = truth_llh
    out = []
    for t, la, lo, h in rows:
        dn = (la - la0) * 111320
        de = (lo - lo0) * 111320 * np.cos(np.radians(la0))
        du = h - h0
        out.append((t, dn, de, du))
    return out


def hms_to_s(t):
    hms = t.split(' ')[1]
    h, m, s = hms.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)


def summarize(label, d, tail_frac=2 / 3):
    arr = np.array([(x[1], x[2], x[3]) for x in d])
    rms3d = np.sqrt((arr ** 2).sum(axis=1))
    n = len(d)
    tail = rms3d[int(n * (1 - tail_frac)):]
    return (f'{label:28s} n={n:4d}  3D err: initial={rms3d[0]:7.3f}m  '
            f'final={rms3d[-1]:7.3f}m  tail-mean={tail.mean():7.3f}m  tail-max={tail.max():7.3f}m')


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    session_dir = sys.argv[1]
    truth_llh = (float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))
    prefix = sys.argv[5] if len(sys.argv) > 5 else 'gs_static'

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    summary_lines = []

    for label, fn, color in (('rapid precise', f'{prefix}_precise.pos', 'tab:orange'),
                              ('ultra-rapid', f'{prefix}_ultra.pos', 'tab:green')):
        rows = load_pos(f'{session_dir}/{fn}')
        d = diffs(rows, truth_llh)
        line = summarize(label, d)
        print(line)
        summary_lines.append(line)

        t0 = hms_to_s(d[0][0])
        ts = [hms_to_s(x[0]) - t0 for x in d]
        r3 = [np.sqrt(x[1] ** 2 + x[2] ** 2 + x[3] ** 2) for x in d]
        axes[0].plot(ts, r3, '.', ms=3, color=color, label=f'{label} 3D error')
        axes[1].plot(ts, [x[3] for x in d], '.', ms=3, color=color, label=f'{label} up error')

    axes[0].set_ylabel('3D position error vs surveyed truth (m)')
    axes[0].set_title(f'{session_dir}: GS static-PPP position convergence')
    axes[0].legend(fontsize=8)
    axes[0].set_ylim(0, 20)

    axes[1].axhline(0, color='k', lw=0.5)
    axes[1].set_ylabel('Up error (m)')
    axes[1].set_xlabel('Elapsed time (s)')
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    out = f'{session_dir}/{prefix}_compare.png'
    plt.savefig(out, dpi=120)
    print(f'wrote {out}')

    with open(f'{session_dir}/{prefix}_compare.txt', 'w') as f:
        f.write(f'=== {session_dir} static PPP (position estimated) vs surveyed truth {truth_llh} ===\n')
        f.write('\n'.join(summary_lines) + '\n')
    print(f'wrote {session_dir}/{prefix}_compare.txt')


if __name__ == '__main__':
    main()
