#!/usr/bin/env python3
#
#  Pocket SDR - Position Solution Visualizer
#
#  Reads $POS entries from a pocket_trk log file and plots:
#    - Latitude / Longitude / Height time series
#    - Horizontal scatter (deviation from reference)
#    - Solution quality and standard deviation over time
#
#  Usage:
#    python3 pocket_pos_plot.py <logfile> [options]
#
#  Options:
#    --ref LAT,LON,HGT   known reference position (deg, deg, m)
#    --ppp-only          show only PPP (Q=6) epochs
#    --out PNGFILE       save figure instead of displaying
#
import sys, argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

SOLQ_LABEL = {0:'---', 1:'FIX', 2:'FLT', 3:'SBS', 4:'DGP', 5:'SPP', 6:'PPP', 7:'DR'}
SOLQ_COLOR = {0:'gray', 1:'blue', 2:'cyan', 3:'magenta', 4:'orange', 5:'green', 6:'red', 7:'purple'}

def parse_log(path):
    """Parse $POS lines; return structured array."""
    rows = []
    with open(path) as f:
        for line in f:
            if not line.startswith('$POS,'):
                continue
            parts = line.strip().split(',')
            if len(parts) < 17:
                continue
            try:
                t    = float(parts[1])
                lat  = float(parts[8])
                lon  = float(parts[9])
                hgt  = float(parts[10])
                q    = int(parts[11])
                ns   = int(parts[12])
                stdn = float(parts[13])
                stde = float(parts[14])
                stdu = float(parts[15])
                dtr  = float(parts[16])
                rows.append((t, lat, lon, hgt, q, ns, stdn, stde, stdu, dtr))
            except (ValueError, IndexError):
                continue
    dtype = [('t','f8'),('lat','f8'),('lon','f8'),('hgt','f8'),
             ('q','i4'),('ns','i4'),('stdn','f8'),('stde','f8'),
             ('stdu','f8'),('dtr','f8')]
    return np.array(rows, dtype=dtype)

def latlon_to_ne(lat, lon, lat0, lon0):
    """Approximate N/E offset in metres from reference (lat0, lon0) in deg."""
    R = 6378137.0
    dlat = np.deg2rad(lat - lat0)
    dlon = np.deg2rad(lon - lon0)
    N = dlat * R
    E = dlon * R * np.cos(np.deg2rad(lat0))
    return N, E

def main():
    ap = argparse.ArgumentParser(description='Pocket SDR position log visualizer')
    ap.add_argument('logfile', help='pocket_trk log file')
    ap.add_argument('--ref', default=None,
                    help='reference position LAT,LON,HGT (deg,deg,m)')
    ap.add_argument('--ppp-only', action='store_true',
                    help='show only PPP (Q=6) epochs')
    ap.add_argument('--out', default=None, help='save figure to PNG file')
    args = ap.parse_args()

    data = parse_log(args.logfile)
    if len(data) == 0:
        print('No $POS records found in', args.logfile)
        sys.exit(1)

    if args.ppp_only:
        data = data[data['q'] == 6]
        if len(data) == 0:
            print('No PPP (Q=6) epochs found.')
            sys.exit(1)

    # reference position
    if args.ref:
        try:
            ref = [float(x) for x in args.ref.split(',')]
            lat0, lon0, hgt0 = ref
        except Exception:
            print('--ref format: LAT,LON,HGT'); sys.exit(1)
    else:
        # use mean of all valid epochs as reference
        lat0 = np.mean(data['lat'])
        lon0 = np.mean(data['lon'])
        hgt0 = np.mean(data['hgt'])
        print(f'Reference (mean): {lat0:.9f}  {lon0:.9f}  {hgt0:.3f}')

    t = (data['t'] - data['t'][0]) / 60.0  # minutes from start
    N, E = latlon_to_ne(data['lat'], data['lon'], lat0, lon0)
    U = data['hgt'] - hgt0

    q_vals = np.unique(data['q'])

    # ── layout ─────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f'PocketSDR Position Solution — {args.logfile}', fontsize=12)
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.42, wspace=0.35)

    ax_n   = fig.add_subplot(gs[0, 0])
    ax_e   = fig.add_subplot(gs[1, 0])
    ax_u   = fig.add_subplot(gs[2, 0])
    ax_ne  = fig.add_subplot(gs[0:2, 1], aspect='equal')
    ax_std = fig.add_subplot(gs[2, 1])

    def scatter_by_q(ax, x, y, **kw):
        for q in q_vals:
            m = data['q'] == q
            ax.scatter(x[m], y[m], s=4, c=SOLQ_COLOR.get(q,'gray'),
                       label=SOLQ_LABEL.get(q, str(q)), **kw)

    # N, E, U time series
    for ax, series, ylabel in [
            (ax_n, N, 'North (m)'),
            (ax_e, E, 'East (m)'),
            (ax_u, U, 'Up (m)')]:
        scatter_by_q(ax, t, series)
        ax.axhline(0, color='k', lw=0.5, ls='--')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Time (min)')
        ax.grid(True, ls=':', lw=0.5)
        rms = np.sqrt(np.mean(series**2))
        ax.set_title(f'RMS {rms:.3f} m', fontsize=9)

    # Horizontal scatter
    scatter_by_q(ax_ne, E, N)
    ax_ne.axhline(0, color='k', lw=0.5, ls='--')
    ax_ne.axvline(0, color='k', lw=0.5, ls='--')
    ax_ne.set_xlabel('East (m)')
    ax_ne.set_ylabel('North (m)')
    ax_ne.set_title('Horizontal scatter')
    ax_ne.grid(True, ls=':', lw=0.5)
    hrms = np.sqrt(np.mean(N**2 + E**2))
    ax_ne.set_title(f'Horizontal scatter  (HRMS={hrms:.3f} m)')

    # Standard deviation over time
    ax_std.plot(t, data['stdn'], label='σN', lw=0.8)
    ax_std.plot(t, data['stde'], label='σE', lw=0.8)
    ax_std.plot(t, data['stdu'], label='σU', lw=0.8)
    ax_std.set_xlabel('Time (min)')
    ax_std.set_ylabel('Std dev (m)')
    ax_std.set_title('Formal accuracy (σ)')
    ax_std.legend(fontsize=8, loc='upper right')
    ax_std.set_ylim(bottom=0)
    ax_std.grid(True, ls=':', lw=0.5)

    # Legend for solution quality
    legend_handles = [Line2D([0],[0], marker='o', color='w',
                             markerfacecolor=SOLQ_COLOR.get(q,'gray'),
                             markersize=6, label=SOLQ_LABEL.get(q,str(q)))
                      for q in q_vals]
    fig.legend(handles=legend_handles, title='Quality', loc='lower center',
               ncol=len(q_vals), fontsize=9, frameon=True,
               bbox_to_anchor=(0.5, 0.01))

    # Summary text
    total = len(data)
    ppp_n = np.sum(data['q'] == 6)
    spp_n = np.sum(data['q'] == 5)
    print(f'Epochs: {total}  SPP: {spp_n}  PPP: {ppp_n}  ({100*ppp_n/total:.1f}% PPP)')
    print(f'N  RMS={np.sqrt(np.mean(N**2)):.3f} m  '
          f'E  RMS={np.sqrt(np.mean(E**2)):.3f} m  '
          f'U  RMS={np.sqrt(np.mean(U**2)):.3f} m')
    if ppp_n > 0:
        ppp = data[data['q']==6]
        Np, Ep = latlon_to_ne(ppp['lat'], ppp['lon'], lat0, lon0)
        Up = ppp['hgt'] - hgt0
        print(f'PPP only — N RMS={np.sqrt(np.mean(Np**2)):.3f} m  '
              f'E RMS={np.sqrt(np.mean(Ep**2)):.3f} m  '
              f'U RMS={np.sqrt(np.mean(Up**2)):.3f} m')

    if args.out:
        plt.savefig(args.out, dpi=150, bbox_inches='tight')
        print('saved:', args.out)
    else:
        plt.show()

if __name__ == '__main__':
    main()
