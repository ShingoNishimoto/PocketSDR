#!/usr/bin/env python3
#
#  pocket_rtk_postproc.py - One-shot RTK post-processing wrapper.
#
#  Wraps the manual pipeline built up over the drone-experiment debugging
#  session into a single command:
#    1. rover.obs  <- pocket_obs2rnx.py pocket.log   (offline, no internet)
#    2. rover.nav  <- pocket_eph2rnx.py pocket.log   (offline, no internet --
#       broadcast ephemeris is identical everywhere on Earth, so there is no
#       need to wait on GA's nav product, which can lag a day or more behind
#       the flight; this sidesteps that wait entirely)
#    3. base obs   <- fetched from the GA GNSS Data Centre API (needs internet;
#       this is the one thing that's genuinely specific to the base station
#       and can't be generated locally) -- skipped if already present.
#    4. rnx2rtkp   <- run with conf/rtk_kinematic.conf
#
#  Usage:
#    python3 pocket_rtk_postproc.py work/20260705162229
#    python3 pocket_rtk_postproc.py work/20260705162229 --skip-fetch
#    python3 pocket_rtk_postproc.py work/20260705162229 \
#        --rnx2rtkp ../RTKLIB_fork/app/consapp/rnx2rtkp/gcc/rnx2rtkp
#
#  IMPORTANT: the Ubuntu/Debian `rtklib` apt package (2.4.3.b34+dfsg-1) does
#  not correctly honor `ant2-postype=rinexhead` -- it silently leaves the base
#  station position at (6378137,0,0), which fails the elevation mask for
#  every satellite and produces a header-only .pos file with zero epochs, no
#  error reported. You need a build from https://github.com/rtklibexplorer/RTKLIB
#  (or current github.com/tomojitakasu/RTKLIB) -- see app/consapp/rnx2rtkp.
#  This script warns if it falls back to whatever "rnx2rtkp" resolves to on
#  PATH, since that is very likely the broken apt build.
#
import sys, os, re, argparse, subprocess, glob
from datetime import datetime, timedelta, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pocket_ga_rinex_fetch as fetch

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CONF = os.path.join(REPO_ROOT, 'conf', 'rtk_kinematic.conf')
# Common location if you cloned demo5 as a sibling of this repo, per the
# session's own setup: ../RTKLIB_fork/app/consapp/rnx2rtkp/gcc/rnx2rtkp
DEFAULT_RNX2RTKP_CANDIDATES = [
    os.path.join(os.path.dirname(REPO_ROOT), 'RTKLIB_fork', 'app', 'consapp',
                 'rnx2rtkp', 'gcc', 'rnx2rtkp'),
]

# GA marks base obs files with "_MO." (rinex3 "M: mixed" Observation) in the
# name; rover.obs (from pocket_obs2rnx.py) never matches this, so a glob on
# it is a safe way to find already-downloaded base chunks without needing to
# track exact filenames returned by a previous fetch.
BASE_OBS_GLOB = '*_MO.*'


def run(cmd, **kw):
    print(f"$ {' '.join(cmd)}")
    return subprocess.run(cmd, **kw)


def find_rnx2rtkp(explicit):
    if explicit:
        if not os.path.isfile(explicit):
            print(f'error: --rnx2rtkp {explicit} not found', file=sys.stderr)
            sys.exit(1)
        return explicit
    for cand in DEFAULT_RNX2RTKP_CANDIDATES:
        if os.path.isfile(cand):
            return cand
    import shutil
    on_path = shutil.which('rnx2rtkp')
    if on_path:
        print(f'WARNING: no demo5 build found; falling back to {on_path!r} on '
              'PATH. If that resolves to the Ubuntu/Debian apt package '
              '(2.4.3.b34+dfsg-1), ant2-postype=rinexhead will NOT work and '
              'you will silently get zero solution epochs. Build from '
              'https://github.com/rtklibexplorer/RTKLIB (app/consapp/rnx2rtkp) '
              'and pass --rnx2rtkp <path> instead.', file=sys.stderr)
        return on_path
    print('error: no rnx2rtkp found (checked demo5 default location and '
          'PATH). Pass --rnx2rtkp <path>.', file=sys.stderr)
    sys.exit(1)


def existing_base_obs(session_dir):
    return sorted(p for p in glob.glob(os.path.join(session_dir, BASE_OBS_GLOB))
                  if not p.endswith(('.obs', '.nav')))


# GA filenames encode their own coverage window, e.g.
# STR100AUS_S_20261890015_15M_01S_MO.crx -> year 2026, day-of-year 189,
# 00:15 start, 15-minute period. Used to detect when *some* base obs file
# exists but doesn't actually span the rover's full flight window (e.g. only
# one of two needed 15-min chunks got fetched) -- silently processing with a
# partial base file produces a large stretch of unreferenced rover epochs,
# not an error, so this has to be checked explicitly rather than assumed.
_GA_FNAME_RE = re.compile(r'_(\d{4})(\d{3})(\d{4})_(\d{2})([DHM])_')
_PERIOD_DUR = {('01', 'D'): timedelta(days=1), ('01', 'H'): timedelta(hours=1),
               ('15', 'M'): timedelta(minutes=15)}


def ga_file_coverage(path):
    """GA-style RINEX filename -> (start, end) UTC coverage, or None if the
    name doesn't match (e.g. a manually-supplied file) -- callers treat None
    as "assume it covers whatever's needed" rather than fail closed."""
    m = _GA_FNAME_RE.search(os.path.basename(path))
    if not m:
        return None
    year, doy, hhmm, num, unit = m.groups()
    dur = _PERIOD_DUR.get((num, unit))
    if dur is None:
        return None
    start = (datetime(int(year), 1, 1, tzinfo=timezone.utc) +
             timedelta(days=int(doy) - 1, hours=int(hhmm[:2]), minutes=int(hhmm[2:])))
    return start, start + dur


def covers_range(paths, start, end):
    """True if paths' merged GA-filename coverage spans [start, end].

    Files with unparseable names are ignored for the check (assumed fine);
    if none of the files parse, the check can't be done at all and this
    conservatively returns False so a fetch is attempted.
    """
    windows = sorted(w for w in (ga_file_coverage(p) for p in paths) if w)
    if not windows:
        return False
    merged = [windows[0]]
    for s, e in windows[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return any(s <= start and e >= end for s, e in merged)


def main():
    ap = argparse.ArgumentParser(
        description='Regenerate rover.obs/rover.nav from pocket.log, fetch '
                    'base obs if needed, and run rnx2rtkp -- one command.')
    ap.add_argument('session_dir', help='directory containing pocket.log')
    ap.add_argument('--conf', default=DEFAULT_CONF,
                    help=f'rnx2rtkp config file (default: {DEFAULT_CONF})')
    ap.add_argument('--rnx2rtkp', default=None,
                    help='path to rnx2rtkp binary (default: auto-detect, see '
                         'module docstring for why this matters)')
    ap.add_argument('--station', default='STR1',
                    help='base station id (default: STR1, Mount Stromlo/Canberra)')
    ap.add_argument('--fallback-station', default='TID1',
                    help='fallback base station id (default: TID1, Tidbinbilla)')
    ap.add_argument('--period', default='15M', choices=['01D', '01H', '15M'],
                    help='base obs file period (default: 15M, high-rate)')
    ap.add_argument('--skip-fetch', action='store_true',
                    help='use already-downloaded base obs in session_dir; '
                         "don't hit the network at all")
    ap.add_argument('--out', default='rtk.pos',
                    help='output .pos filename, relative to session_dir '
                         '(default: rtk.pos)')
    args = ap.parse_args()

    session = args.session_dir.rstrip('/')
    log_path = os.path.join(session, 'pocket.log')
    if not os.path.isfile(log_path):
        print(f'error: {log_path} not found', file=sys.stderr)
        sys.exit(1)

    rover_obs = os.path.join(session, 'rover.obs')
    rover_nav = os.path.join(session, 'rover.nav')

    print('== 1/4: rover.obs (offline) ==')
    run([sys.executable, os.path.join(REPO_ROOT, 'python', 'pocket_obs2rnx.py'),
        log_path, '-o', rover_obs], check=True)

    print('\n== 2/4: rover.nav (offline) ==')
    run([sys.executable, os.path.join(REPO_ROOT, 'python', 'pocket_eph2rnx.py'),
        log_path, '-o', rover_nav], check=False)
    if not os.path.isfile(rover_nav):
        print(f'WARNING: {rover_nav} was not produced -- likely no $EPH '
              'records in this log (older flights predate the log_mask fix '
              'in src/sdr_func.c). Falling back to a base-station nav file '
              'if one is already present in the session directory.',
              file=sys.stderr)

    print('\n== 3/4: base observations ==')
    base_obs = existing_base_obs(session)
    if not args.skip_fetch:
        start, end = fetch.flight_range_from_log(log_path)
    if base_obs and not args.skip_fetch and covers_range(base_obs, start, end):
        print(f'Existing base obs ({len(base_obs)} file(s)) already cover the '
              'full flight window; skipping fetch.')
    elif args.skip_fetch:
        if base_obs:
            print(f'Using {len(base_obs)} existing base obs file(s) as-is '
                  '(--skip-fetch); coverage not checked.')
    else:
        if base_obs:
            print(f'Existing base obs ({len(base_obs)} file(s)) do not fully '
                  f'cover {start.isoformat()}..{end.isoformat()} -- fetching '
                  'to fill the gap (existing files are kept, not overwritten).')
        obs_records = fetch.query_api([args.station], start, end, args.period,
                                       ['obs'], '3')
        if not obs_records and args.fallback_station:
            print(f'No obs for {args.station}; trying {args.fallback_station}')
            obs_records = fetch.query_api([args.fallback_station], start, end,
                                           args.period, ['obs'], '3')
        if not obs_records:
            print('error: no base obs files found for this window',
                  file=sys.stderr)
            sys.exit(1)
        for rec in obs_records:
            path = fetch.download(rec, session)
            path = fetch.maybe_gunzip(path)
            fetch.maybe_crx2rnx(path)
        base_obs = existing_base_obs(session)
    if not base_obs:
        print('error: no base obs files available (pass --skip-fetch only if '
              'you already have some, or check network access)', file=sys.stderr)
        sys.exit(1)
    print(f'Base obs: {len(base_obs)} file(s)')

    # A single glob argument lets rnx2rtkp itself expand+concatenate
    # multi-chunk base data; quoting matters so the shell doesn't pre-expand it.
    base_arg = os.path.join(session, BASE_OBS_GLOB) if len(base_obs) > 1 \
        else base_obs[0]

    nav_arg = rover_nav if os.path.isfile(rover_nav) else None
    if nav_arg is None:
        base_nav = sorted(glob.glob(os.path.join(session, '*_MN.*')))
        if not base_nav:
            print('error: no nav source available (neither rover.nav nor a '
                  'downloaded base _MN nav file)', file=sys.stderr)
            sys.exit(1)
        nav_arg = base_nav[0]
    print(f'Nav source: {nav_arg}')

    print('\n== 4/4: rnx2rtkp ==')
    rnx2rtkp = find_rnx2rtkp(args.rnx2rtkp)
    out_path = os.path.join(session, args.out)
    run([rnx2rtkp, '-k', args.conf, rover_obs, base_arg, nav_arg,
        '-o', out_path], check=True)

    with open(out_path) as f:
        rows = [l for l in f if not l.startswith('%')]
    print(f'\n{out_path}: {len(rows)} solution epoch(s)')
    if rows:
        from collections import Counter
        q = Counter(l.split()[5] for l in rows if len(l.split()) > 5)
        for k in sorted(q):
            print(f'  Q={k}: {q[k]}')
        print(f'\nNext: python3 python/pocket_pos_plot.py {out_path} '
              '--ahd LAT LON AHD N --ant 0.5')
    else:
        print('WARNING: zero solution epochs -- see the module docstring for '
              'the ant2-postype/rnx2rtkp-build gotcha.', file=sys.stderr)


if __name__ == '__main__':
    main()
