#!/usr/bin/env python3
"""Skyplot (azimuth/elevation track) and CN0-vs-time plot (L1 + L2) for all
satellites tracked in a pocket.log, from its "$SAT" and "$OBS" lines:

    $SAT,time,year,month,day,hour,min,sec,sat,pvt,obs,cn0,az,el,res
        time  receiver time (s elapsed)
        sat   satellite ID (e.g. G07, J02)
        pvt   PVT status (0: not used, 1: used in solution)
        obs   L1 obs data status (0: not available, 1: available)
        cn0   L1 signal C/N0 (dB-Hz)
        az,el azimuth/elevation (deg)
        res   L1 pseudorange residual (m)
(see out_log_sat() in src/sdr_pvt.c -- used for the skyplot, az/el only)

    $OBS,time,year,month,day,hour,min,sec,sat,code,cn0,pr,cp,dop,lli,fcn,ch
        code  RINEX obs code (e.g. "1C" for L1CA, "2S" for L2CM)
        cn0   C/N0 for that signal (dB-Hz)
(see out_log_obs() in src/sdr_pvt.c -- used for the CN0 plot, both
L1 (1C) and L2 (2S) per satellite, since $SAT only carries L1)

Usage:
    python3 pocket_sat_plot.py <pocket.log> [out_prefix]

out_prefix defaults to the log's own directory + "sat"; writes
<out_prefix>_skyplot.png and <out_prefix>_cn0.png.
"""
import sys
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SAT_RE = re.compile(
    r'\$SAT,([\d.]+),(\d+),(\d+),(\d+),(\d+),(\d+),([\d.]+),(\w+),(\d),(\d),'
    r'([\d.\-]+),([\d.\-]+),([\d.\-]+),([\d.\-]+)')
OBS_RE = re.compile(
    r'\$OBS,([\d.]+),(\d+),(\d+),(\d+),(\d+),(\d+),([\d.]+),(\w+),(\w+),'
    r'([\d.\-]+),')

# RINEX obs code -> frequency label, for the CN0 plot
CODE_FREQ = {'1C': 'L1', '2S': 'L2', '2C': 'L2', '2L': 'L2', '2X': 'L2'}


def parse_sat_lines(log_path):
    """Returns {sat: {'t': [...], 'az': [...], 'el': [...], 'cn0': [...],
    'pvt': [...]}} sorted by time within each satellite."""
    data = {}
    with open(log_path, errors='replace') as f:
        for line in f:
            m = SAT_RE.search(line)
            if not m:
                continue
            (t, y, mo, d, h, mi, s, sat, pvt, obs, cn0, az, el, res) = m.groups()
            rec = data.setdefault(sat, {'t': [], 'az': [], 'el': [], 'cn0': [], 'pvt': []})
            rec['t'].append(float(t))
            rec['az'].append(float(az))
            rec['el'].append(float(el))
            rec['cn0'].append(float(cn0))
            rec['pvt'].append(int(pvt))
    for sat, rec in data.items():
        order = np.argsort(rec['t'])
        for k in rec:
            rec[k] = np.array(rec[k])[order]
    return data


def parse_obs_lines(log_path):
    """Returns {sat: {'L1': {'t':[...], 'cn0':[...]}, 'L2': {...}}}, sorted
    by time within each satellite/frequency. Codes not in CODE_FREQ (e.g.
    L5) are ignored -- this script only plots L1/L2 per the request."""
    data = {}
    with open(log_path, errors='replace') as f:
        for line in f:
            m = OBS_RE.search(line)
            if not m:
                continue
            t, y, mo, d, h, mi, s, sat, code, cn0 = m.groups()
            freq = CODE_FREQ.get(code)
            if freq is None:
                continue
            rec = data.setdefault(sat, {}).setdefault(freq, {'t': [], 'cn0': []})
            rec['t'].append(float(t))
            rec['cn0'].append(float(cn0))
    for sat, freqs in data.items():
        for freq, rec in freqs.items():
            order = np.argsort(rec['t'])
            rec['t'] = np.array(rec['t'])[order]
            rec['cn0'] = np.array(rec['cn0'])[order]
    return data


def sat_colors(sats):
    cmap = plt.get_cmap('tab20')
    return {sat: cmap(i % 20) for i, sat in enumerate(sorted(sats))}


def plot_skyplot(data, colors, out_path, title):
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlim(90, 0)  # elevation: 90 at center, 0 at edge
    ax.set_rticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_rlabel_position(135)

    for sat in sorted(data):
        rec = data[sat]
        el_mask = rec['el'] > 0
        if not el_mask.any():
            continue
        theta = np.radians(rec['az'][el_mask])
        r = rec['el'][el_mask]
        ax.plot(theta, r, '-', lw=1.2, color=colors[sat], label=sat)
        ax.plot(theta[-1], r[-1], 'o', ms=5, color=colors[sat])
        ax.annotate(sat, (theta[-1], r[-1]), fontsize=8, color=colors[sat],
                    xytext=(3, 3), textcoords='offset points')

    ax.set_title(title)
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path)
    print(f'wrote {out_path}')


FREQ_MARKER = {'L1': '.', 'L2': 'x'}
FREQ_MS = {'L1': 2, 'L2': 3}


def plot_cn0(obs_data, colors, out_path, title):
    fig, ax = plt.subplots(figsize=(12, 7))
    t0 = min(rec['t'].min() for freqs in obs_data.values() for rec in freqs.values())
    for sat in sorted(obs_data):
        for freq in ('L1', 'L2'):
            rec = obs_data[sat].get(freq)
            if rec is None:
                continue
            ax.plot(rec['t'] - t0, rec['cn0'], FREQ_MARKER[freq], ms=FREQ_MS[freq],
                    color=colors[sat], label=sat if freq == 'L1' else None)
    ax.set_xlabel('Elapsed time (s)')
    ax.set_ylabel('C/N0 (dB-Hz)')
    ax.set_title(title)
    sat_legend = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=8, ncol=1,
                            title='satellite')
    ax.add_artist(sat_legend)
    freq_handles = [Line2D([0], [0], marker=FREQ_MARKER[f], color='k', ls='',
                            ms=FREQ_MS[f], label=f) for f in ('L1', 'L2')]
    ax.legend(handles=freq_handles, loc='lower left', bbox_to_anchor=(1.01, 0.0),
              fontsize=8, title='freq')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path)
    print(f'wrote {out_path}')


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    log_path = sys.argv[1]
    if len(sys.argv) > 2:
        out_prefix = sys.argv[2]
    else:
        out_prefix = log_path.rsplit('/', 1)[0] + '/sat' if '/' in log_path else 'sat'

    data = parse_sat_lines(log_path)
    obs_data = parse_obs_lines(log_path)
    if not data:
        print(f'no $SAT lines found in {log_path}')
        sys.exit(1)
    if not obs_data:
        print(f'no $OBS lines found in {log_path}')
        sys.exit(1)
    print(f'{log_path}: {len(data)} satellites tracked: {sorted(data)}')
    n_l2 = sum(1 for freqs in obs_data.values() if 'L2' in freqs)
    print(f'{n_l2}/{len(obs_data)} satellites have L2 obs data')

    colors = sat_colors(set(data) | set(obs_data))
    plot_skyplot(data, colors, f'{out_prefix}_skyplot.pdf', 'satellite skyplot')
    plot_cn0(obs_data, colors, f'{out_prefix}_cn0.pdf', 'satellite C/N0 (L1/L2)')


if __name__ == '__main__':
    main()
