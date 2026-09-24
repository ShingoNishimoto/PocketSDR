#!/usr/bin/env python3
"""5-way GS clock comparison: fixed-position vs static (position-estimated)
PPP modes, each with rapid-precise and ultra-rapid ephemeris, plus the
session's own real-time broadcast-ephemeris run. Three rows:
  1. clock offset (ms)
  2. PPP-reported clock std (m)
  3. clock offset difference from a reference series (m) -- default
     reference is fixed-position + rapid-precise, the tightest-converging
     of the five; every series (including the reference, identically zero)
     is interpolated onto the reference's own time grid and differenced in
     metres, since the raw offsets differ only sub-millisecond -- unreadable
     directly against row 1's ms scale.

Prerequisite: run gs_ppp_export.py twice first (--prefix gs and --prefix
gs_static) to produce <session_dir>/gs_ppp_results.xlsx and
gs_static_ppp_results.xlsx (each with RapidPrecise/UltraRapid sheets).

Usage:
    python3 gs_plot_fixed_vs_static.py <session_dir> <brdc_xlsx> [brdc_sheet]

Example:
    python3 gs_plot_fixed_vs_static.py ../work/gs/3-5 pocket_gs.xlsx
    python3 gs_plot_fixed_vs_static.py ../work/gs/3-4 pocket.xlsx
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CLIGHT = 299792458.0
REF_LABEL = 'fixed-pos, rapid precise'

# Sized for a single column of a double-column journal layout (IEEE/Elsevier
# single-column width is ~3.5in/88.9mm) -- figure is built at exactly that
# width so it drops in at 100% scale with no further shrinking (which is
# what usually makes journal figure text illegible).
COL_WIDTH_IN = 3.5
FONT_SIZE = 7
plt.rcParams.update({
    'font.size': FONT_SIZE,
    'axes.titlesize': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE - 1,
    'ytick.labelsize': FONT_SIZE - 1,
    'legend.fontsize': FONT_SIZE - 1,
    'lines.markersize': 2,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
})


def hms_to_s(t):
    h, m, s = t.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)


DTR_OUTLIER_THRESH = 0.0005  # 0.5ms; real epoch-to-epoch drift is ~tens of ns


def clean(df, col):
    d = df[['datetime_gpst', col]].dropna().copy()
    tstr = d['datetime_gpst'].astype(str).str.split('T').str[1]
    d['t_s'] = tstr.apply(hms_to_s)
    d = d.sort_values('t_s')
    if col == 'dtr_s':
        # A handful of rows (mostly at session start, and at data-quality
        # boundaries the pppdbg merge doesn't cover -- see gs_ppp_export.py)
        # fall back to parse_rtklib_pos_as_ref()'s crude ~1ms timestamp-
        # fraction approximation of dtr_s instead of the real PPP-estimated
        # value, showing up as multi-millisecond spikes/drops against the
        # otherwise smooth ~tens-of-microseconds drift trend. Drop anything
        # further than DTR_OUTLIER_THRESH from the series median.
        med = d[col].median()
        d = d[(d[col] - med).abs() <= DTR_OUTLIER_THRESH]
    return d


def tail_stats(label, dtr, std):
    n = len(dtr)
    dtr_tail = dtr.iloc[n // 3:]
    std_tail = std.iloc[len(std) // 3:]
    print(f'{label:34s} n={n:5d}  dtr std(tail)={dtr_tail["dtr_s"].std()*1e3:8.4f} ms  '
          f'mean clk_std(tail)={std_tail["ppp_clk_std_m"].mean():7.4f} m')


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    session_dir, brdc_xlsx = sys.argv[1], sys.argv[2]
    brdc_sheet = sys.argv[3] if len(sys.argv) > 3 else 'PPP'

    fixed_precise = pd.read_excel(f'{session_dir}/gs_ppp_results.xlsx', sheet_name='RapidPrecise')
    fixed_ultra = pd.read_excel(f'{session_dir}/gs_ppp_results.xlsx', sheet_name='UltraRapid')
    static_precise = pd.read_excel(f'{session_dir}/gs_static_ppp_results.xlsx', sheet_name='RapidPrecise')
    static_ultra = pd.read_excel(f'{session_dir}/gs_static_ppp_results.xlsx', sheet_name='UltraRapid')
    brdc = pd.read_excel(f'{session_dir}/{brdc_xlsx}', sheet_name=brdc_sheet)

    series = {
        'broadcast (real-time)': (brdc, 'tab:blue'),
        REF_LABEL: (fixed_precise, 'tab:orange'),
        'fixed-pos, ultra-rapid': (fixed_ultra, 'tab:green'),
        'static (pos est), rapid precise': (static_precise, 'tab:red'),
        'static (pos est), ultra-rapid': (static_ultra, 'tab:purple'),
    }
    # Short legend text -- full names are long enough that even a 2-column
    # legend wraps/collides at 3.5in; full names still go in the printed
    # tail-stats table and belong in the figure caption.
    DISPLAY = {
        'broadcast (real-time)': 'Fixed+Broadcast',
        REF_LABEL: 'Fixed+Rapid (ref.)',
        'fixed-pos, ultra-rapid': 'Fixed+Ultra',
        'static (pos est), rapid precise': 'Static+Rapid',
        'static (pos est), ultra-rapid': 'Static+Ultra',
    }
    cleaned = {k: (clean(df, 'dtr_s'), clean(df, 'ppp_clk_std_m'), color)
               for k, (df, color) in series.items()}
    t0 = min(dtr['t_s'].min() for dtr, _, _ in cleaned.values())

    # Clip every series to the shortest common window (see gs_ppp_export.py's
    # module docstring for why the precise/ultra-rapid runs can stop early).
    t_end = min(dtr['t_s'].max() for dtr, _, _ in cleaned.values())
    print(f'clipping comparison plot to common window: t0={t0:.1f} t_end={t_end:.1f} '
          f'({t_end - t0:.1f}s, vs broadcast\'s own '
          f'{cleaned["broadcast (real-time)"][0]["t_s"].max() - t0:.1f}s)')
    cleaned = {k: (dtr[dtr['t_s'] <= t_end], std[std['t_s'] <= t_end], color)
               for k, (dtr, std, color) in cleaned.items()}

    ref_dtr, _, _ = cleaned[REF_LABEL]
    ref_t = ref_dtr['t_s'].to_numpy()
    ref_y = ref_dtr['dtr_s'].to_numpy()

    # Tight x-range: the actual plotted data (post outlier-filter, post
    # common-window clip), not matplotlib's default 5%-padded autoscale.
    t_lo = min(min(dtr['t_s'].min(), std['t_s'].min()) for dtr, std, _ in cleaned.values()) - t0
    t_hi = max(max(dtr['t_s'].max(), std['t_s'].max()) for dtr, std, _ in cleaned.values()) - t0
    # t_hi = 1160 # for 3-5

    with plt.rc_context({'font.size': 7, 'axes.labelsize': 5.5, 'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5}):
        fig, axes = plt.subplots(3, 1, figsize=(COL_WIDTH_IN, 3.5), sharex=True,
                                gridspec_kw={'height_ratios': [3, 3, 3], 'hspace': 0.12})
        handles, labels = [], []
        for label, (dtr, std, color) in cleaned.items():
            h, = axes[0].plot(dtr['t_s'] - t0, dtr['dtr_s'] * 1e3, '.', color=color, ms=1)
            handles.append(h)
            labels.append(DISPLAY[label])
            interp_y = np.interp(ref_t, dtr['t_s'].to_numpy(), dtr['dtr_s'].to_numpy())
            if label != REF_LABEL:
                diff_m = (interp_y - ref_y) * CLIGHT
                axes[1].plot(ref_t - t0, diff_m, '.', color=color, ms=1)
            axes[2].plot(std['t_s'] - t0, std['ppp_clk_std_m'], '.', color=color, ms=1)

        axes[0].set_ylabel('Clock offset (ms)')
        axes[1].axhline(0, color='k', lw=0.4)
        axes[1].set_ylabel('Clock diff. vs ref (m)')
        axes[1].set_ylim(-2, 10)
        axes[2].set_ylabel('Clock std (m)')
        axes[2].set_ylim(0, 1.5)
        axes[-1].set_xlabel('Elapsed time (s)')
        for ax in axes:
            ax.set_xlim(t_lo, t_hi)
            ax.margins(x=0)

        # One shared legend for all 3 panels (all use the same 5 series) instead
        # of repeating it per-axes -- saves the vertical space a per-panel
        # legend would eat at this width, which is most of it at 3.5in. Reserve
        # room for it with an explicit subplots_adjust (not tight_layout, which
        # doesn't know how to make space for a figure-level legend and warned/
        # collided with it) sized for 3 legend rows at this font/column count.
        fig.subplots_adjust(left=0.18, right=0.97, bottom=0.10, top=0.98)
        leg = fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.58, 1.0),
                        ncol=2, fontsize=FONT_SIZE - 1, columnspacing=1.0, handletextpad=0.3,
                        labelspacing=0.25, frameon=False)
        # Query the legend's actual rendered height instead of guessing a top
        # margin by hand -- robust to future font/row-count/label-length changes.
        fig.canvas.draw()
        leg_bottom_fig = leg.get_window_extent().transformed(fig.transFigure.inverted()).y0
        fig.subplots_adjust(top=leg_bottom_fig - 0.02)

        out = f'{session_dir}/gs_fixed_vs_static_clock.pdf'
        fig.savefig(out)
        print(f'wrote {out}  (middle-panel reference series: "{REF_LABEL}" -- '
            f'state this in the figure caption, it is no longer spelled out on the axis)')

        print()
        for label, (dtr, std, _) in cleaned.items():
            tail_stats(label, dtr, std)


if __name__ == '__main__':
    main()
