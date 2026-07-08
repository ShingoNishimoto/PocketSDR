#!/usr/bin/env python3
#
#  pocket_pos_plot.py - PocketSDR position accuracy evaluator
#
#  Reads $POS entries from pocket_trk log file(s) OR RTKLIB .pos files,
#  computes N/E/U errors relative to a known reference position, and plots:
#    - Height convergence with +/-1sigma band and reference line
#    - North / East / Up error time series with formal +/-sigma envelope
#    - Horizontal position scatter
#    - Formal standard deviation (sigmaN, sigmaE, sigmaU) over time
#  Prints mean, std, and RMS per component to the terminal.
#
#  File format is auto-detected:
#    pocket.log  — lines starting with $POS (pocket_trk output)
#    *.pos       — RTKLIB rnx2rtkp / rtkpost LLH output (date time lat lon hgt Q ...)
#
#  Usage:
#    python3 pocket_pos_plot.py <file> [file2 ...] [options]
#
#  Reference position options (at least one required for error evaluation):
#    --ref LAT LON HGT   WGS-84 ellipsoidal (deg deg m)
#    --ahd LAT LON AHD N AHD orthometric height + geoid undulation N (m);
#                         WGS-84 hgt = AHD + N + ant.  Obtain N from
#                         Geoscience Australia AUSGeoid2020 online tool.
#  Other options:
#    --ant H              antenna height above ground reference point (m)
#    --ppp-only           include only PPP/FIX (Q=6 for PPP, Q=1 for RTK-FIX)
#    --after T            statistics window: only T seconds after first PPP/FIX
#    --out FILE           save figure to PNG instead of displaying
#
import sys, argparse, os
from datetime import datetime, timezone
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

_WGS84_A  = 6378137.0
_WGS84_E2 = 2 / 298.257223563 - (1 / 298.257223563) ** 2

SOLQ_LABEL  = {0:'---', 1:'FIX', 2:'FLT', 3:'SBS', 4:'DGP',
               5:'SPP', 6:'PPP', 7:'DR'}
SOLQ_COLOR  = {0:'gray',    1:'blue',   2:'cyan',  3:'magenta',
               4:'orange',  5:'green',  6:'red',   7:'purple'}
FILE_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
               'tab:purple', 'tab:brown', 'tab:pink',  'tab:olive']

_DTYPE = [('t','f8'), ('lat','f8'), ('lon','f8'), ('hgt','f8'),
          ('q','i4'), ('ns','i4'),
          ('stdn','f8'), ('stde','f8'), ('stdu','f8')]

# ── parsing ────────────────────────────────────────────────────────────────────

def _parse_pocket_log(path, tag='$POS,'):
    """Parse $POS (or $REFPOS, same field layout) records from a pocket_trk log."""
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            if not line.startswith(tag):
                continue
            p = line.strip().split(',')
            if len(p) < 17:
                continue
            try:
                rows.append((
                    float(p[1]),
                    float(p[8]), float(p[9]), float(p[10]),
                    int(p[11]),  int(p[12]),
                    float(p[13]), float(p[14]), float(p[15]),
                ))
            except (ValueError, IndexError):
                continue
    return np.array(rows, dtype=_DTYPE)


def _gpst_to_ts(date_str, time_str):
    """Convert RTKLIB GPST 'YYYY/MM/DD HH:MM:SS.sss' to float Unix timestamp."""
    dt = datetime.strptime(f'{date_str} {time_str}', '%Y/%m/%d %H:%M:%S.%f')
    return dt.replace(tzinfo=timezone.utc).timestamp()


def _parse_rtklib_pos(path):
    """Parse RTKLIB rnx2rtkp/rtkpost LLH .pos file.

    Expected (space-separated, lines not starting with %):
      date  time  lat  lon  hgt  Q  ns  sdn  sde  sdu  ...
    """
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%') or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 10:
                continue
            try:
                t  = _gpst_to_ts(p[0], p[1])
                rows.append((
                    t,
                    float(p[2]), float(p[3]), float(p[4]),
                    int(p[5]),   int(p[6]),
                    float(p[7]), float(p[8]), float(p[9]),
                ))
            except (ValueError, IndexError):
                continue
    return np.array(rows, dtype=_DTYPE)


def _detect_format(path):
    """Return 'pocket' or 'rtklib' based on file content."""
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('$POS,'):
                return 'pocket'
            if line and not line.startswith('%') and not line.startswith('#'):
                p = line.split()
                if len(p) >= 10 and '/' in p[0] and ':' in p[1]:
                    return 'rtklib'
    return 'pocket'


def parse_log(path, tag='$POS,'):
    """Return structured numpy array — auto-detects pocket.log vs RTKLIB .pos.

    tag selects '$POS,' (main solver, default) or '$REFPOS,' (reference
    solver) when path is a pocket.log; ignored for RTKLIB .pos files.
    """
    fmt = _detect_format(path)
    if fmt == 'rtklib':
        data = _parse_rtklib_pos(path)
        if len(data) > 0:
            # Normalise t to seconds-from-start (same as pocket.log field p[1])
            data['t'] = data['t'] - data['t'][0]
        return data
    return _parse_pocket_log(path, tag)

# ── geodesy ───────────────────────────────────────────────────────────────────

def llh_to_ecef(lat_deg, lon_deg, hgt_m):
    """WGS-84 geodetic (deg, deg, m) → ECEF (m), vectorised."""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    N = _WGS84_A / np.sqrt(1.0 - _WGS84_E2 * np.sin(lat)**2)
    X = (N + hgt_m) * np.cos(lat) * np.cos(lon)
    Y = (N + hgt_m) * np.cos(lat) * np.sin(lon)
    Z = (N * (1.0 - _WGS84_E2) + hgt_m) * np.sin(lat)
    return X, Y, Z

def ecef_to_neu(dx, dy, dz, lat0_deg, lon0_deg):
    """ECEF delta → North/East/Up at reference (lat0, lon0), vectorised."""
    phi = np.deg2rad(lat0_deg)
    lam = np.deg2rad(lon0_deg)
    sp, cp = np.sin(phi), np.cos(phi)
    sl, cl = np.sin(lam), np.cos(lam)
    n = -sp*cl*dx - sp*sl*dy + cp*dz
    e = -sl*dx       + cl*dy
    u =  cp*cl*dx + cp*sl*dy + sp*dz
    return n, e, u

def compute_neu(data, lat0, lon0, hgt0):
    """Compute N/E/U errors (m) for all rows in data vs reference."""
    X0, Y0, Z0 = llh_to_ecef(lat0, lon0, hgt0)
    X, Y, Z = llh_to_ecef(data['lat'], data['lon'], data['hgt'])
    return ecef_to_neu(X - X0, Y - Y0, Z - Z0, lat0, lon0)

# ── statistics ────────────────────────────────────────────────────────────────

def print_stats(label, Ne, Ee, Ue):
    H2 = Ne**2 + Ee**2
    P2 = H2 + Ue**2
    hdr = f'{"":16s} {"N (m)":>9s} {"E (m)":>9s} {"U (m)":>9s}' \
          f' {"H2D (m)":>9s} {"3D (m)":>9s}'
    mean_row = (f'{"Mean (bias)":16s} {np.mean(Ne):+9.3f} {np.mean(Ee):+9.3f}'
                f' {np.mean(Ue):+9.3f}')
    std_row  = (f'{"Std (1sigma)":16s} {np.std(Ne):9.3f}  {np.std(Ee):9.3f}'
                f'  {np.std(Ue):9.3f}')
    rms_row  = (f'{"RMS":16s} {np.sqrt(np.mean(Ne**2)):9.3f}'
                f'  {np.sqrt(np.mean(Ee**2)):9.3f}'
                f'  {np.sqrt(np.mean(Ue**2)):9.3f}'
                f'  {np.sqrt(np.mean(H2)):9.3f}'
                f'  {np.sqrt(np.mean(P2)):9.3f}')
    print(f'\n{label}  (n = {len(Ne)})')
    print(hdr)
    print(mean_row)
    print(std_row)
    print(rms_row)

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description='PocketSDR position accuracy evaluator')
    ap.add_argument('logfiles', nargs='+', help='pocket_trk log file(s)')
    ap.add_argument('--ref', nargs=3, type=float,
                    metavar=('LAT', 'LON', 'HGT'),
                    help='WGS-84 reference: LAT LON HGT (deg deg m ellipsoidal)')
    ap.add_argument('--ahd', nargs=4, type=float,
                    metavar=('LAT', 'LON', 'AHD', 'N'),
                    help='AHD reference: LAT LON AHD N  (WGS-84 = AHD+N+ant)')
    ap.add_argument('--ant', type=float, default=0.0,
                    help='antenna height above ground reference point (m, default 0)')
    ap.add_argument('--ppp-only', action='store_true',
                    help='include only best-quality epochs: Q=6 (PPP) or Q=1 (RTK-FIX)')
    ap.add_argument('--after', type=float, default=0.0,
                    help='statistics window: seconds after first PPP epoch (default=all)')
    ap.add_argument('--refpos', action='store_true',
                    help='plot $REFPOS (reference solver) instead of $POS '
                         '(main solver); ignored for RTKLIB .pos files')
    ap.add_argument('--out', default=None,
                    help='save figure to file (.png or .pdf, by extension)')
    args = ap.parse_args()
    tag = '$REFPOS,' if args.refpos else '$POS,'

    # ── resolve reference position ────────────────────────────────────────────
    lat0 = lon0 = hgt0 = None

    if args.ref:
        lat0, lon0, hgt0 = args.ref
        if args.ant != 0.0:
            print(f'Note: --ant {args.ant} m ignored for --ref (HGT should '
                  'already be at antenna phase centre)')

    elif args.ahd:
        lat0, lon0, ahd, N_geoid = args.ahd
        hgt0 = ahd + N_geoid + args.ant
        print(f'Reference: {lat0:.7f}  {lon0:.7f}')
        print(f'  AHD {ahd:.3f} m + geoid N {N_geoid:.3f} m'
              f' + ant {args.ant:.3f} m = {hgt0:.3f} m WGS-84 ellipsoidal')

    have_ref = (lat0 is not None)

    # ── load data ─────────────────────────────────────────────────────────────
    datasets = []
    for path in args.logfiles:
        d = parse_log(path, tag)
        if len(d) == 0:
            print(f'Warning: no {tag.rstrip(",")} records in {path}')
            continue
        if args.ppp_only:
            d = d[(d['q'] == 6) | (d['q'] == 1)]  # PPP or RTK-FIX
            if len(d) == 0:
                print(f'Warning: no PPP/FIX epochs in {path}')
                continue
        datasets.append((os.path.basename(path), d))

    if not datasets:
        print('No valid data.')
        sys.exit(1)

    if not have_ref:
        d0 = datasets[0][1]
        ppp = d0[d0['q'] == 6]
        ref_data = ppp if len(ppp) > 0 else d0
        lat0 = float(np.mean(ref_data['lat']))
        lon0 = float(np.mean(ref_data['lon']))
        hgt0 = float(np.mean(ref_data['hgt']))
        print('No reference given — using mean of first file:')
        print(f'  lat={lat0:.9f}  lon={lon0:.9f}  hgt={hgt0:.3f} m WGS-84')

    # ── figure layout ─────────────────────────────────────────────────────────
    multi = len(datasets) > 1
    title = ', '.join(n for n, _ in datasets)
    fig = plt.figure(figsize=(14, 12))
    fig.suptitle(f'PocketSDR Position Accuracy — {title}', fontsize=11)

    gs = gridspec.GridSpec(4, 2, figure=fig,
                           hspace=0.50, wspace=0.32,
                           height_ratios=[1.5, 1, 1, 1])
    ax_h   = fig.add_subplot(gs[0, 0])         # Height convergence (absolute)
    ax_n   = fig.add_subplot(gs[1, 0])         # North error
    ax_e   = fig.add_subplot(gs[2, 0])         # East error
    ax_u   = fig.add_subplot(gs[3, 0])         # Up error
    ax_ne  = fig.add_subplot(gs[0:2, 1],       # Horizontal scatter
                             aspect='equal')
    ax_std = fig.add_subplot(gs[2:4, 1])       # Formal sigma time series

    # ── per-dataset plotting ──────────────────────────────────────────────────
    for idx, (name, data) in enumerate(datasets):
        fcolor = FILE_COLORS[idx % len(FILE_COLORS)]
        t = (data['t'] - data['t'][0]) / 60.0      # minutes from session start
        Ne, Ee, Ue = compute_neu(data, lat0, lon0, hgt0)
        q_vals = np.unique(data['q'])

        # Height convergence: absolute height with +/-1 sigmaU shading
        ax_h.fill_between(t, data['hgt'] - data['stdu'],
                             data['hgt'] + data['stdu'],
                          alpha=0.2, color=fcolor)
        for q in q_vals:
            m = data['q'] == q
            c   = fcolor if multi else SOLQ_COLOR.get(q, 'gray')
            lbl = (f'{name} ' if multi else '') + SOLQ_LABEL.get(q, str(q))
            ax_h.scatter(t[m], data['hgt'][m], s=3, c=c, label=lbl, zorder=3)

        # N / E / U error time series with formal sigma shading
        for ax, err, std_key in [(ax_n, Ne, 'stdn'),
                                  (ax_e, Ee, 'stde'),
                                  (ax_u, Ue, 'stdu')]:
            sig = data[std_key]
            ax.fill_between(t, -sig, +sig, alpha=0.15, color='gray')
            for q in q_vals:
                m = data['q'] == q
                c = fcolor if multi else SOLQ_COLOR.get(q, 'gray')
                ax.scatter(t[m], err[m], s=3, c=c, zorder=3)

        # Horizontal scatter
        for q in q_vals:
            m = data['q'] == q
            c   = fcolor if multi else SOLQ_COLOR.get(q, 'gray')
            lbl = (f'{name} ' if multi else '') + SOLQ_LABEL.get(q, str(q))
            ax_ne.scatter(Ee[m], Ne[m], s=5, c=c, label=lbl, zorder=3)

        # Formal sigma time series
        ls = ['-', '--', ':', '-.'][idx % 4]
        pfx = f'{name} ' if multi else ''
        ax_std.plot(t, data['stdn'], ls=ls, color='C0', lw=0.9,
                    label=f'{pfx}sigmaN')
        ax_std.plot(t, data['stde'], ls=ls, color='C1', lw=0.9,
                    label=f'{pfx}sigmaE')
        ax_std.plot(t, data['stdu'], ls=ls, color='C2', lw=0.9,
                    label=f'{pfx}sigmaU')

        # Statistics: best-quality epochs (PPP Q=6, or RTK-FIX Q=1)
        fix_mask = (data['q'] == 6) | (data['q'] == 1)   # PPP or RTK-FIX
        ppp_mask = data['q'] == 6                          # PPP only (for label)
        n_total  = len(data)
        n_fix    = int(np.sum(fix_mask))
        n_spp    = int(np.sum(data['q'] == 5))
        fix_label = 'PPP' if np.any(ppp_mask) else 'FIX'
        print(f'\n{name}: {n_total} epochs  '
              f'SPP={n_spp}  {fix_label}={n_fix} ({100*n_fix/n_total:.1f}%)')

        if np.any(fix_mask):
            t0_ppp   = data['t'][np.argmax(fix_mask)]
            stat_mask = fix_mask
            stat_label = f'{name} — {fix_label} all'
            if args.after > 0:
                stat_mask = fix_mask & (data['t'] >= t0_ppp + args.after)
                stat_label = f'{name} — {fix_label} after {args.after:.0f} s'
                # mark the cutoff on time-series panels
                t_cut = (t0_ppp + args.after - data['t'][0]) / 60.0
                for ax in (ax_h, ax_n, ax_e, ax_u):
                    ax.axvline(t_cut, color='k', lw=0.7, ls=':', alpha=0.5)
            if np.sum(stat_mask) > 0:
                print_stats(stat_label,
                            Ne[stat_mask], Ee[stat_mask], Ue[stat_mask])
            else:
                print(f'  (no epochs in stats window — '
                      'try reducing --after or check PPP coverage)')
        else:
            print_stats(f'{name} — all epochs', Ne, Ee, Ue)

    # ── axes decoration ───────────────────────────────────────────────────────
    if have_ref:
        ax_h.axhline(hgt0, color='k', lw=1.2, ls='--',
                     label=f'reference {hgt0:.2f} m', zorder=4)
    ax_h.set_ylabel('Height WGS-84 (m)')
    ax_h.set_xlabel('Time (min)')
    ax_h.set_title('Height convergence  (shading = ±1σU)')
    ax_h.grid(True, ls=':', lw=0.5)
    ax_h.legend(fontsize=7, loc='upper right')

    for ax, ylabel, std_sym in [(ax_n, 'North error (m)', 'σN'),
                                 (ax_e, 'East error (m)',  'σE'),
                                 (ax_u, 'Up error (m)',    'σU')]:
        ax.axhline(0, color='k', lw=0.8, ls='--')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Time (min)')
        ax.grid(True, ls=':', lw=0.5)
        ax.text(0.01, 0.97, f'shading = ±{std_sym}',
                transform=ax.transAxes, fontsize=7, va='top', color='gray')

    ax_ne.axhline(0, color='k', lw=0.5, ls='--')
    ax_ne.axvline(0, color='k', lw=0.5, ls='--')
    ax_ne.set_xlabel('East error (m)')
    ax_ne.set_ylabel('North error (m)')
    ax_ne.set_title('Horizontal scatter')
    ax_ne.grid(True, ls=':', lw=0.5)

    ax_std.set_xlabel('Time (min)')
    ax_std.set_ylabel('Formal std dev (m)')
    ax_std.set_title('Formal accuracy  (σN, σE, σU)')
    ax_std.set_ylim(bottom=0)
    ax_std.legend(fontsize=7, loc='upper right')
    ax_std.grid(True, ls=':', lw=0.5)

    # Solution quality legend (single-file only; multi-file uses ax_h legend)
    if not multi:
        _, data = datasets[0]
        q_vals = np.unique(data['q'])
        handles = [Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=SOLQ_COLOR.get(q, 'gray'),
                          markersize=6, label=SOLQ_LABEL.get(q, str(q)))
                   for q in q_vals]
        fig.legend(handles=handles, title='Solution Q', loc='lower center',
                   ncol=len(q_vals), fontsize=9, frameon=True,
                   bbox_to_anchor=(0.5, 0.005))
    else:
        ax_ne.legend(fontsize=7, loc='upper right')

    if args.out:
        plt.savefig(args.out, dpi=150, bbox_inches='tight')
        print(f'\nFigure saved: {args.out}')
    else:
        plt.show()

if __name__ == '__main__':
    main()
