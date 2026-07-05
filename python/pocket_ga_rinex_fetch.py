#!/usr/bin/env python3
#
#  pocket_ga_rinex_fetch.py - Fetch base-station RINEX obs/nav files from the
#  Geoscience Australia GNSS Data Centre API for RTK post-processing.
#
#  API reference: https://data.gnss.ga.gov.au/docs/rinex-file-query/v1.0/web-api-access.html
#    GET /api/rinexFiles?stationId=...&startDate=...&endDate=...&filePeriod=...
#                        &fileType=...&rinexVersion=...
#  Response is a JSON array of file records, each with a "fileLocation" S3/HTTPS
#  download URL. No authentication needed for open (unrestricted) stations.
#
#  Usage:
#    # derive the exact flight window automatically from a pocket.log session
#    # (fetches only the 15-min high-rate chunks the flight actually spans):
#    python3 pocket_ga_rinex_fetch.py --from-log work/20260705.../pocket.log \
#        --out-dir work/20260705...
#
#    # or give the whole UTC day explicitly (grabs every 15-min chunk that day):
#    python3 pocket_ga_rinex_fetch.py --date 2026-07-05 --station STR1 \
#        --fallback-station TID1 --out-dir base_rinex
#
#  Note: GA only offers a nav product at 01D (and 01H) periods, never at 15M
#  (confirmed live: querying filePeriod=15M&fileType=nav 404s) -- nav is
#  always fetched separately at 01D regardless of --period, since one day's
#  broadcast nav file covers the whole day's flight windows anyway.
#
#  Then, if the flight spans more than one 15-min chunk, pass all the base
#  obs chunks as a single quoted wildcard argument -- rnx2rtkp expands and
#  concatenates same-station multi-file obs input internally:
#    rnx2rtkp -k rtk_kinematic.conf rover.obs "<out-dir>/STR1*_MO.crx.gz" \
#        <out-dir>/STR1*_MN.rnx.gz -o rtk.pos
#
import sys, os, argparse, json, gzip, shutil, subprocess, urllib.request, urllib.parse
from datetime import datetime, timedelta, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pocket_export as pe

API_BASE = 'https://data.gnss.ga.gov.au/api/rinexFiles'
NAV_PERIOD = '01D'  # GA only serves nav products at 01D/01H, never 15M


def iso(dt):
    return dt.strftime('%Y-%m-%dT%H:%M:%SZ')


def row_to_dt(r):
    """$POS row's GPST calendar fields -> UTC-labelled datetime (no leap shift,
    consistent with how the GA API timestamps its own file windows)."""
    si = int(r['sec'])
    extra_min, si = divmod(si, 60)
    return (datetime(r['year'], r['month'], r['day'], r['hour'], r['min'], 0,
                     tzinfo=timezone.utc) + timedelta(minutes=extra_min, seconds=si))


def flight_range_from_log(path):
    """Exact flight start/end from the first and last $POS epoch.

    The GA API matches by containment (verified live: a query range only
    needs to touch a file's window to include it), so passing the precise
    flight span -- rather than the whole day -- fetches only the 15-min/
    hourly chunks the flight actually needs.
    """
    rows = pe.parse_log(path)
    if not rows:
        print(f'No $POS records found in {path}', file=sys.stderr)
        sys.exit(1)
    return row_to_dt(rows[0]), row_to_dt(rows[-1])


def whole_day_range(year, month, day):
    """(00:00:00, 23:59:59) UTC -- matches every chunk of any period that day."""
    start = datetime(year, month, day, 0, 0, 0, tzinfo=timezone.utc)
    end = datetime(year, month, day, 23, 59, 59, tzinfo=timezone.utc)
    return start, end


def query_api(station_ids, start, end, period, filetypes, rinex_version):
    params = {
        'stationId': ','.join(station_ids),
        'startDate': iso(start),
        'endDate': iso(end),
        'filePeriod': period,
        'fileType': ','.join(filetypes),
        'rinexVersion': str(rinex_version),
    }
    url = API_BASE + '?' + urllib.parse.urlencode(params)
    print(f'Query: {url}')
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            body = resp.read()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # The API 404s (rather than 200 + []) when nothing matches --
            # normal for very recent dates where daily/hourly products
            # (nav in particular) haven't been assembled yet.
            body = e.read()
        else:
            raise
    records = json.loads(body)
    # Some records may be metadata-only error entries (no fileLocation)
    return [r for r in records if r.get('fileLocation')]


def download(record, out_dir):
    url = record['fileLocation']
    name = os.path.basename(urllib.parse.urlparse(url).path)
    dest = os.path.join(out_dir, name)
    print(f'  {record.get("siteId","?"):5s} {record.get("fileType","?"):4s} '
          f'{record.get("filePeriod","?"):4s} rinex{record.get("rinexVersion","?")} '
          f'{record.get("fileSize", 0):>10d} bytes -> {dest}')
    urllib.request.urlretrieve(url, dest)
    return dest


def maybe_gunzip(path):
    if not path.endswith('.gz'):
        return path
    out_path = path[:-3]
    with gzip.open(path, 'rb') as fin, open(out_path, 'wb') as fout:
        shutil.copyfileobj(fin, fout)
    return out_path


def is_hatanaka(path):
    """RINEX2 compact obs: '.YYd'/'.YYD' (e.g. str10010.19d). RINEX3: '.crx'."""
    base = os.path.basename(path)
    if base.lower().endswith('.crx'):
        return True
    _, ext = os.path.splitext(base)
    return len(ext) == 4 and ext[1:3].isdigit() and ext[3] in ('d', 'D')


def maybe_crx2rnx(path):
    if not is_hatanaka(path):
        return path
    # RTKLIB's own rtk_uncompress() shells out to lowercase 'crx2rnx' on Unix;
    # check both spellings since the official binary is distributed either way.
    exe = shutil.which('crx2rnx') or shutil.which('CRX2RNX')
    if exe is None:
        print(f'  NOTE: {path} is Hatanaka-compressed; no crx2rnx/CRX2RNX on '
              'PATH. Either install it, or -- since this repo\'s RTKLIB '
              'shells out to "crx2rnx" internally -- just point rnx2rtkp at '
              'the original .crx.gz file directly and it will decompress '
              'on the fly, provided gzip and crx2rnx are on PATH.')
        return path
    subprocess.run([exe, '-f', path], check=True)
    if path.lower().endswith('.crx'):
        out_path = path[:-4] + '.rnx'
    else:
        out_path = path[:-1] + ('O' if path[-1] == 'D' else 'o')
    return out_path if os.path.exists(out_path) else path


def main():
    ap = argparse.ArgumentParser(
        description='Fetch base-station RINEX obs/nav from the GA GNSS Data Centre API')
    src = ap.add_mutually_exclusive_group(required=False)
    src.add_argument('--from-log', metavar='pocket.log',
                     help='derive the exact flight window to fetch from a pocket.log '
                          'session (recommended -- fetches only the needed chunks)')
    src.add_argument('--date', metavar='YYYY-MM-DD',
                     help='whole UTC calendar day to fetch (every chunk that day)')
    ap.add_argument('--start', help='explicit ISO8601 start, e.g. 2026-07-05T00:00:00Z '
                    '(overrides --date/--from-log)')
    ap.add_argument('--end', help='explicit ISO8601 end (paired with --start)')
    ap.add_argument('--station', default='STR1',
                    help='comma-separated station id(s), tried together (default STR1, '
                         'Mount Stromlo/Canberra)')
    ap.add_argument('--fallback-station', default='TID1',
                    help='station id to retry with if --station returns nothing '
                         '(default TID1, Tidbinbilla)')
    ap.add_argument('--period', default='15M', choices=['01D', '01H', '15M'],
                    help='obs file period (default 15M, high-rate 1Hz for the '
                         'kinematic rover baseline; kept online 1 year. '
                         '01H is 30s data kept 14 days; 01D is 30s, kept forever. '
                         'Nav is always fetched at 01D regardless of this setting.')
    ap.add_argument('--filetype', default='obs,nav',
                    help='comma-separated file types (default obs,nav)')
    ap.add_argument('--rinex-version', default='3',
                    help='RINEX version filter (default 3; this repo\'s RTKLIB '
                         'fork does not parse RINEX 4)')
    ap.add_argument('--out-dir', default='.', help='download directory (default .)')
    ap.add_argument('--no-decompress', action='store_true',
                    help='skip gunzip / Hatanaka decompression of downloaded files')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.start:
        start = datetime.strptime(args.start, '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
        end = datetime.strptime(args.end, '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
    elif args.from_log:
        start, end = flight_range_from_log(args.from_log)
    elif args.date:
        y, m, d = (int(x) for x in args.date.split('-'))
        start, end = whole_day_range(y, m, d)
    else:
        ap.error('one of --from-log, --date, or --start/--end is required')

    print(f'Fetch window (UTC): {iso(start)} .. {iso(end)}')

    filetypes = set(args.filetype.split(','))
    stations = args.station.split(',')

    def fetch(types, period):
        if not types:
            return []
        recs = query_api(stations, start, end, period, types, args.rinex_version)
        if not recs and args.fallback_station:
            print(f'No {period} {sorted(types)} files for station(s) {stations}; '
                  f'retrying with fallback {args.fallback_station!r}')
            recs = query_api([args.fallback_station], start, end, period,
                             types, args.rinex_version)
        return recs

    # nav has no 15M product on GA -- always fetch it at NAV_PERIOD separately
    obs_types = filetypes - {'nav'}
    obs_records = fetch(obs_types, args.period)
    nav_records = fetch(filetypes & {'nav'}, NAV_PERIOD)
    records = obs_records + nav_records

    if not records:
        print('No matching RINEX files found.', file=sys.stderr)
        sys.exit(1)

    if 'nav' in filetypes and not nav_records:
        print(f'\nWARNING: no {NAV_PERIOD} nav file found for this window at any '
              f'station tried. GA daily/hourly products are usually only '
              f'assembled after the UTC day (or hour) has fully elapsed -- '
              f'if {iso(end)} is very recent, this is most likely data '
              f'latency, not a real gap. Try again later, or fetch nav for '
              f'an earlier date if you just need something to test against.',
              file=sys.stderr)

    print(f'\nDownloading {len(records)} file(s) to {args.out_dir}/')
    obs_paths, nav_paths = [], []
    for rec in records:
        path = download(rec, args.out_dir)
        if not args.no_decompress:
            path = maybe_gunzip(path)
            path = maybe_crx2rnx(path)
        (obs_paths if rec.get('fileType') == 'obs' else nav_paths).append(path)

    print('\nDone:')
    for p in obs_paths + nav_paths:
        print(f'  {p}')

    print('\nNext step:')
    if len(obs_paths) > 1:
        # Multiple time-chunked base obs files (e.g. several 15-min high-rate
        # files spanning a flight >15 min): pass them as ONE quoted wildcard
        # argument -- rnx2rtkp expands and concatenates same-station obs
        # input internally, rather than passing each file as a separate arg.
        prefix = os.path.commonprefix([os.path.basename(p) for p in obs_paths])
        obs_arg = f'"{args.out_dir}/{prefix}*"' if prefix else \
            '"' + '" "'.join(obs_paths) + '"'
        print(f'  {len(obs_paths)} base obs chunks -- pass as one wildcard arg:')
    else:
        obs_arg = obs_paths[0] if obs_paths else '<base_obs>'
    print(f'  rnx2rtkp -k rtk_kinematic.conf rover.obs {obs_arg} '
          f'{nav_paths[0] if nav_paths else "<base_nav>"} -o rtk.pos')


if __name__ == '__main__':
    main()
