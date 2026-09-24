#!/usr/bin/env python3
#
#  pocket_export.py - Export PocketSDR pocket.log to geo/data formats.
#
#  Reads $POS (and optionally $REFPOS) entries and writes any combination of:
#    KML     (Google Earth, PPP and SPP as separate folders)
#    GPX     (GPS track format, PPP and SPP as separate tracks)
#    GeoJSON (web mapping, per-point properties)
#    CSV     (spreadsheet-compatible)
#    Excel   (.xlsx with All / PPP / SPP sheets; Ref / RefPPP / RefSPP added when
#             $REFPOS records are present)
#
#  Usage:
#    python3 pocket_export.py pocket.log [--all] [--kml F] [--gpx F]
#                             [--geojson F] [--csv F] [--excel F]
#                             [--ref-kml F] [--ref-gpx F] [--ref-geojson F]
#                             [--ref-csv F] [--ref-excel F]
#
#  With --all, all formats are written; file names are derived from the log
#  file base name (e.g. pocket.log → pocket.kml, pocket.gpx, ...).
#  Reference solution files are written automatically when $REFPOS records
#  exist (pocket_ref.kml, pocket_ref.csv, …).
#
import sys, argparse, json, os
from datetime import datetime, timezone, timedelta
import numpy as np
import pandas as pd

# WGS84 ellipsoid parameters
_WGS84_A  = 6378137.0
_WGS84_E2 = 2 / 298.257223563 - (1 / 298.257223563) ** 2

SOLQ_LABEL = {0:'NONE', 1:'FIX', 2:'FLOAT', 3:'SBAS', 4:'DGPS',
              5:'SPP', 6:'PPP', 7:'DR'}

# GPS is ahead of UTC by 18 leap seconds as of 2017-01-01; update if needed.
GPS_LEAP_SECONDS = 18

# ── parsing ────────────────────────────────────────────────────────────────────

def parse_log(path):
    """Return list of dicts parsed from $POS lines in a pocket.log file.

    $POS format:
      t,year,month,day,hour,min,sec,lat,lon,hgt,q,ns,stdn,stde,stdu,dtr [17 fields]
      + ecef_x,ecef_y,ecef_z,vel_x,vel_y,vel_z                          [+6, len>=23]
      + gdop,pdop,hdop,vdop                                              [+4, len>=27]
      + dtr_drift                                                        [+1, len>=28]
    """
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$POS,'):
                continue
            p = line.split(',')
            if len(p) < 17:
                continue
            try:
                r = {
                    't':     float(p[1]),
                    'year':  int(float(p[2])), 'month': int(float(p[3])),
                    'day':   int(float(p[4])), 'hour':  int(float(p[5])),
                    'min':   int(float(p[6])), 'sec':   float(p[7]),
                    'lat':   float(p[8]),  'lon':  float(p[9]),
                    'hgt':   float(p[10]),
                    'q':     int(p[11]),   'ns':   int(p[12]),
                    'stdn':  float(p[13]), 'stde': float(p[14]),
                    'stdu':  float(p[15]), 'dtr':  float(p[16]),
                }
                if len(p) >= 23:
                    r['ecef_x'] = float(p[17])
                    r['ecef_y'] = float(p[18])
                    r['ecef_z'] = float(p[19])
                    r['vel_x']  = float(p[20])
                    r['vel_y']  = float(p[21])
                    r['vel_z']  = float(p[22])
                else:
                    r['ecef_x'], r['ecef_y'], r['ecef_z'] = \
                        llh_to_ecef(r['lat'], r['lon'], r['hgt'])
                    r['vel_x'] = r['vel_y'] = r['vel_z'] = float('nan')
                if len(p) >= 27:
                    r['gdop'] = float(p[23])
                    r['pdop'] = float(p[24])
                    r['hdop'] = float(p[25])
                    r['vdop'] = float(p[26])
                if len(p) >= 28:
                    r['dtr_drift'] = float(p[27])
                # lsq() in RTKLIB's pntpos() can leave NaN in a failed solve's
                # sol.dtr[]/rr[] while still reporting stat!=0 (seen live: a
                # singular normal-equations matrix from too few sats). Such
                # rows aren't usable position fixes -- skip them here so every
                # consumer (calendar conversion, plotting, comparison) doesn't
                # have to separately defend against NaN calendar fields.
                if not all(np.isfinite([r['sec'], r['lat'], r['lon'], r['hgt'], r['dtr']])):
                    continue
                rows.append(r)
            except (ValueError, IndexError):
                continue
    return rows


def parse_sc_aowr_log(path):
    """Parse $LOG,<t>,SC_AOWR entries and return a dict keyed by t.

    Line format:  $LOG,<t>,SC_AOWR clk=<float> drift=<float> dt=<float> hist=<int>
    Returns: {t: {'clk': float, 'drift': float, 'dt': float, 'hist': int}, ...}
    """
    result = {}
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$LOG,'):
                continue
            p = line.split(',', 2)
            if len(p) < 3 or not p[2].startswith('SC_AOWR '):
                continue
            try:
                t     = float(p[1])
                msg   = p[2]
                clk   = float(msg.split('clk=')[1].split()[0])
                drift = float(msg.split('drift=')[1].split()[0])
                dt    = float(msg.split('dt=')[1].split()[0])
                hist  = int(msg.split('hist=')[1].split()[0])
                result[t] = {'clk': clk, 'drift': drift, 'dt': dt, 'hist': hist}
            except (ValueError, IndexError):
                continue
    return result


def parse_pppos_dbg_log(path):
    """Parse $LOG,<t>,PPPOS_DBG entries for PPP convergence monitoring.

    Handles both old format (no clk_std/zwd/amb_std fields) and new format.
    Returns: {t: {pos_std, stat, ndual, clk, clk_std, zwd, zwd_std,
                  amb_std, nconv, res_c}, ...}
    """
    def _fv(msg, key):
        try:
            return float(msg.split(key + '=')[1].split()[0])
        except (IndexError, ValueError):
            return float('nan')

    def _iv(msg, key):
        try:
            return int(float(msg.split(key + '=')[1].split()[0]))
        except (IndexError, ValueError):
            return None

    result = {}
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$LOG,'):
                continue
            p = line.split(',', 2)
            if len(p) < 3 or not p[2].startswith('PPPOS_DBG '):
                continue
            try:
                t   = float(p[1])
                msg = p[2]
                result[t] = {
                    'pos_std': _fv(msg, 'pos_std'),
                    'stat':    _iv(msg, 'stat'),
                    'ndual':   _iv(msg, 'ndual'),
                    'clk':     _fv(msg, 'clk'),
                    'clk_std': _fv(msg, 'clk_std'),
                    'zwd':     _fv(msg, 'zwd'),
                    'zwd_std': _fv(msg, 'zwd_std'),
                    'amb_std': _fv(msg, 'amb_std'),
                    'nconv':   _iv(msg, 'nconv'),
                    'res_c':   _fv(msg, 'res_c'),
                }
            except (ValueError, IndexError):
                continue
    return result


def parse_spp_seed_log(path, tag):
    """Parse $LOG entries with a given tag and return a dict keyed by t.

    Handles tags emitted by sdr_pvt.c:
      SPP_SEED     – main SPP seed (AOWR-corrected obs; dtr=0 when clock fixed)
      REF_SPP_SEED – reference SPP seed in PPP mode (uncorrected obs, free clock)
      REF_SPP      – reference SPP in SPP mode (uncorrected obs, free clock)

    Line format:  $LOG,<t>,<tag> stat=<int> dtr=<float> [msg=<str>]
    Returns: {t: {'stat': int, 'dtr': float}, ...}
    """
    result = {}
    prefix = tag + ' '
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$LOG,'):
                continue
            p = line.split(',', 2)
            if len(p) < 3 or not p[2].startswith(prefix):
                continue
            try:
                t    = float(p[1])
                msg  = p[2]
                stat = int(msg.split('stat=')[1].split()[0])
                dtr  = float(msg.split('dtr=')[1].split()[0])
                result[t] = {'stat': stat, 'dtr': dtr}
            except (ValueError, IndexError):
                continue
    return result


def parse_refpos_log(path):
    """Return list of dicts parsed from $REFPOS lines (reference/uncorrected solver).

    $REFPOS format (4-unknown solve, clock estimated):
      t,year,month,day,hour,min,sec,lat,lon,hgt,stat,ns,stdn,stde,stdu,dtr [17 fields]
      + gdop,pdop,hdop,vdop                                                 [+4, len>=21]
      + dtr_drift                                                           [+1, len>=22]
      + vel_x,vel_y,vel_z (ECEF m/s)                                       [+3, len>=25]
    ECEF position is computed from lat/lon/hgt (not in the log format).
    """
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$REFPOS,'):
                continue
            p = line.split(',')
            if len(p) < 17:
                continue
            try:
                r = {
                    't':     float(p[1]),
                    'year':  int(float(p[2])), 'month': int(float(p[3])),
                    'day':   int(float(p[4])), 'hour':  int(float(p[5])),
                    'min':   int(float(p[6])), 'sec':   float(p[7]),
                    'lat':   float(p[8]),  'lon':  float(p[9]),
                    'hgt':   float(p[10]),
                    'q':     int(p[11]),   'ns':   int(p[12]),
                    'stdn':  float(p[13]), 'stde': float(p[14]),
                    'stdu':  float(p[15]), 'dtr':  float(p[16]),
                }
                r['ecef_x'], r['ecef_y'], r['ecef_z'] = \
                    llh_to_ecef(r['lat'], r['lon'], r['hgt'])
                r['vel_x'] = r['vel_y'] = r['vel_z'] = float('nan')
                if len(p) >= 21:
                    r['gdop'] = float(p[17])
                    r['pdop'] = float(p[18])
                    r['hdop'] = float(p[19])
                    r['vdop'] = float(p[20])
                if len(p) >= 22:
                    r['dtr_drift'] = float(p[21])
                if len(p) >= 25:
                    r['vel_x'] = float(p[22])
                    r['vel_y'] = float(p[23])
                    r['vel_z'] = float(p[24])
                # See the matching skip in parse_log(): a failed pntpos() lsq()
                # solve can leave NaN in dtr/rr while stat!=0.
                if not all(np.isfinite([r['sec'], r['lat'], r['lon'], r['hgt'], r['dtr']])):
                    continue
                rows.append(r)
            except (ValueError, IndexError):
                continue
    return rows


def parse_dtr_debug_log(path):
    """Parse a DTRDBG debug trace (temporary rtkpos.c instrumentation -- see
    parse_rtklib_pos_as_ref()) into {datetime: dtr} for the forward pass only.

    Line format: DTRDBG,<gps week>,<gps tow>,<sol.stat>,<sol.dtr[0] in seconds>
    printed to stderr right after pppos() returns, once per epoch per pass.
    For pos1-soltype=combined (forward+backward), the trace contains the
    forward pass first (tow increasing) followed by the backward pass (tow
    decreasing) -- only the forward-pass run is kept here, because combres()
    in postpos.c builds the combined output's sols.dtr[0] from solf[i].dtr[0]
    (the forward pass value is carried through unchanged; only position/
    covariance get smoothed) whenever forward and backward cover the same
    epoch, which is the common case for a continuous observation file.
    """
    GPS_EPOCH = datetime(1980, 1, 6)
    out = {}
    last_tow = -1e18
    with open(path, errors='replace') as f:
        for raw in f:
            for line in raw.split('\r'):
                line = line.strip()
                if not line.startswith('DTRDBG,'):
                    continue
                p = line.split(',')
                if len(p) != 5:
                    continue
                try:
                    week, tow, dtr = int(p[1]), float(p[2]), float(p[4])
                except ValueError:
                    continue
                if tow <= last_tow:
                    break  # backward pass started
                last_tow = tow
                dt = GPS_EPOCH + timedelta(seconds=week * 604800 + tow)
                out[dt] = dtr
    return out


def parse_rtklib_pos_as_ref(path, dtr_debug_path=None):
    """Return list of dicts, in the same shape as parse_refpos_log(), parsed
    from an rnx2rtkp RTKLIB .pos file (llh format, out-timeform=hms/out-
    timesys=gpst -- see conf/rtk_*.conf).

    Used to substitute an offline-reprocessed solution (e.g. standalone PPP
    post-processing, pos1-posmode=ppp-kine) for the onboard reference
    solver's own $REFPOS data, when the latter is unusable (e.g. stuck in
    SPP for most of a session -- see AOWR PPP KF reset bug in sdr_pvt.c).

    dtr_debug_path: optional DTRDBG trace (see parse_dtr_debug_log()) giving
    the true per-epoch receiver clock bias straight from the KF state,
    matched to each .pos row by exact GPST timestamp. Without it, dtr is only
    approximated from the .pos timestamp's own clock-bias-driven sub-second
    offset (see below) -- that approximation is accurate to just ~1 ms (the
    .pos file's own timestamp precision, out-timendec=3), nowhere near the
    ns-level agreement expected between two independent solvers estimating
    the same physical receiver clock. Verified on 2026-07-15 against a
    session with a healthy onboard $REFPOS: DTRDBG-sourced dtr agreed with
    the onboard main solver's independent dtr_s to within ~1e-7 s at the
    same epoch, vs. ~1e-3 s off (and briefly wrong-signed) for the
    timestamp-only approximation.

    Caveats vs a genuine $REFPOS row:
      - Without dtr_debug_path, dtr is reconstructed from the .pos
        timestamp: pntpos() sets sol->time = timeadd(obs[0].time,
        -x[3]/CLIGHT) (note the MINUS -- pntpos.c:433/464), and pppos()
        never touches sol.time afterward (rtkpos.c calls pntpos() to seed
        rtk->sol before pppos(), and only pntpos assigns sol.time), so the
        printed timestamp is raw_obs_epoch MINUS the SPP-seed clock bias,
        i.e. dtr = raw_obs_epoch - reported_time (opposite sign from the
        onboard $REFPOS timestamp convention, time2epoch(timeadd(rsol->time,
        rdtr), ep) in sdr_pvt.c, which adds dtr). PocketSDR rover.obs is
        always sampled on an exact 1 Hz grid, so dtr is recovered from the
        reported timestamp's sub-second offset from the nearest whole
        second -- but see the ms-level precision caveat above; prefer
        dtr_debug_path when available.
      - 't' (elapsed seconds) is relative to this file's own first epoch,
        not pocket.log's own elapsed-time clock, since a standalone rnx2rtkp
        run has no relationship to the onboard process's own start time.
        Calendar fields (year/month/day/hour/min/sec) are the real GPST time
        and unaffected by this.
    """
    _nan = float('nan')
    dtr_debug = parse_dtr_debug_log(dtr_debug_path) if dtr_debug_path else None
    n_debug_miss = 0
    rows = []
    t0 = None
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            p = line.split()
            if len(p) < 10:
                continue
            try:
                dt = datetime.strptime(f'{p[0]} {p[1]}', '%Y/%m/%d %H:%M:%S.%f')
                if t0 is None:
                    t0 = dt
                if dtr_debug is not None:
                    if dt in dtr_debug:
                        dtr = dtr_debug[dt]
                    else:
                        n_debug_miss += 1
                        frac = dt.microsecond / 1e6
                        dtr = -frac if frac <= 0.5 else 1.0 - frac
                else:
                    frac = dt.microsecond / 1e6
                    dtr = -frac if frac <= 0.5 else 1.0 - frac
                r = {
                    't':     (dt - t0).total_seconds(),
                    'year':  dt.year, 'month': dt.month, 'day': dt.day,
                    'hour':  dt.hour, 'min':   dt.minute,
                    'sec':   dt.second + dt.microsecond / 1e6,
                    'lat':   float(p[2]), 'lon': float(p[3]), 'hgt': float(p[4]),
                    'q':     int(p[5]),   'ns':  int(p[6]),
                    'stdn':  float(p[7]), 'stde': float(p[8]), 'stdu': float(p[9]),
                    'dtr':   dtr,  # DTRDBG-sourced if available; else timestamp-approximated
                }
                r['ecef_x'], r['ecef_y'], r['ecef_z'] = \
                    llh_to_ecef(r['lat'], r['lon'], r['hgt'])
                r['vel_x'] = r['vel_y'] = r['vel_z'] = _nan
                if not all(np.isfinite([r['sec'], r['lat'], r['lon'], r['hgt']])):
                    continue
                rows.append(r)
            except (ValueError, IndexError):
                continue
    if dtr_debug is not None and n_debug_miss:
        print(f'Warning: {n_debug_miss}/{len(rows)} rows in {path} had no exact '
              f'DTRDBG timestamp match (fell back to timestamp-approximated dtr '
              '-- likely backward-pass-only epochs near a data gap)',
              file=sys.stderr)
    return rows


def gpst_to_utc(r):
    """Convert $POS/$REFPOS epoch fields (GPST) to UTC datetime object."""
    si = int(r['sec'])
    us = int((r['sec'] - si) * 1e6)
    # time2epoch() can return sec=60.000 at a minute boundary due to float rounding
    extra_min, si = divmod(si, 60)
    dt = datetime(r['year'], r['month'], r['day'],
                  r['hour'], r['min'], si, us, tzinfo=timezone.utc)
    return dt + timedelta(minutes=extra_min) - timedelta(seconds=GPS_LEAP_SECONDS)

def gpst_str(r):
    dt = gpst_to_utc(r) + timedelta(seconds=GPS_LEAP_SECONDS)
    frac = r['sec'] % 1
    return (f"{dt.year}-{dt.month:02d}-{dt.day:02d}"
            f"T{dt.hour:02d}:{dt.minute:02d}:{dt.second + frac:06.3f}")

GPS_EPOCH_TS = datetime(1980, 1, 6, tzinfo=timezone.utc).timestamp()

def gpst_week_tow(r):
    """GPST calendar fields ($POS/$REFPOS, no leap-second shift) -> (week, tow_s)."""
    si = int(r['sec'])
    extra_min, si = divmod(si, 60)
    dt = (datetime(r['year'], r['month'], r['day'], r['hour'], r['min'], 0,
                   tzinfo=timezone.utc) + timedelta(minutes=extra_min, seconds=si))
    total = dt.timestamp() + (r['sec'] - int(r['sec'])) - GPS_EPOCH_TS
    week = int(total // 604800)
    return week, total - week * 604800

def llh_to_ecef(lat_deg, lon_deg, hgt_m):
    """WGS84 geodetic (deg, deg, m) → ECEF (m)."""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    N = _WGS84_A / np.sqrt(1 - _WGS84_E2 * np.sin(lat) ** 2)
    X = (N + hgt_m) * np.cos(lat) * np.cos(lon)
    Y = (N + hgt_m) * np.cos(lat) * np.sin(lon)
    Z = (N * (1 - _WGS84_E2) + hgt_m) * np.sin(lat)
    return X, Y, Z


# ── KML ───────────────────────────────────────────────────────────────────────

# $POS/$REFPOS 'hgt' is WGS-84 *ellipsoidal* height, but Google Earth's
# <altitudeMode>absolute</altitudeMode> is height above mean sea level (the
# EGM96 geoid). The two differ by the local geoid undulation N, which is
# tens of meters and varies slowly with position -- without correcting for
# it, the KML track sits offset vertically from where it should render.
_EGM96_GRID_DIRS = ['/usr/share/proj']


def _egm96_heights(lat_deg, lon_deg, hgt_ellip_m):
    """Convert WGS-84 ellipsoidal heights to EGM96 (MSL) heights via PROJ.

    Returns None (caller falls back to uncorrected ellipsoidal height) if
    pyproj isn't installed or the EGM96 grid can't be found.
    """
    try:
        import pyproj
    except ImportError:
        print('warning: pyproj not installed -- KML altitude will use raw '
              'WGS-84 ellipsoidal height, not EGM96/MSL height. Install '
              "with 'pip install pyproj' to enable the geoid correction.",
              file=sys.stderr)
        return None

    for d in _EGM96_GRID_DIRS:
        if os.path.isdir(d):
            pyproj.datadir.append_data_dir(d)
    pyproj.network.set_network_enabled(False)

    try:
        t = pyproj.Transformer.from_crs('EPSG:4979', 'EPSG:4326+5773',
                                         always_xy=True)
        # sanity check: if the EGM96 grid isn't found, PROJ silently passes
        # the height through unchanged instead of raising.
        _, _, h_test = t.transform(0.0, 0.0, 100.0)
        if abs(h_test - 100.0) < 0.01:
            print('warning: EGM96 geoid grid not found (looked in '
                  f'{_EGM96_GRID_DIRS}) -- KML altitude will use raw '
                  'WGS-84 ellipsoidal height, not EGM96/MSL height.',
                  file=sys.stderr)
            return None
        _, _, hgt_msl = t.transform(lon_deg, lat_deg, hgt_ellip_m)
    except Exception as e:
        print(f'warning: EGM96 geoid transform failed ({e}) -- KML altitude '
              'will use raw WGS-84 ellipsoidal height.', file=sys.stderr)
        return None
    return np.asarray(hgt_msl)


def write_kml(rows, path, label='', use_geoid=True):
    if use_geoid and rows:
        lat = np.array([r['lat'] for r in rows])
        lon = np.array([r['lon'] for r in rows])
        hgt = np.array([r['hgt'] for r in rows])
        hgt_msl = _egm96_heights(lat, lon, hgt)
        if hgt_msl is not None:
            rows = [dict(r, hgt=h) for r, h in zip(rows, hgt_msl)]
            alt_note = 'MSL (EGM96 geoid)'
        else:
            alt_note = 'WGS-84 ellipsoidal (EGM96 geoid correction unavailable)'
    else:
        alt_note = 'WGS-84 ellipsoidal'

    ppp = [r for r in rows if r['q'] == 6]
    spp = [r for r in rows if r['q'] != 6]
    tag = f' ({label})' if label else ''

    def track(pts, name, line_style, dot_style):
        coords = ' '.join(f"{r['lon']},{r['lat']},{r['hgt']}" for r in pts)
        pms = []
        if coords:
            pms.append(
                f'    <Placemark><name>{name} Track</name>'
                f'<styleUrl>#{line_style}</styleUrl>'
                f'<LineString><altitudeMode>absolute</altitudeMode>'
                f'<coordinates>{coords}</coordinates>'
                f'</LineString></Placemark>')
        for r in pts:
            desc = (f"t={r['t']:.1f}s Q={SOLQ_LABEL.get(r['q'],r['q'])} "
                    f"ns={r['ns']} σN={r['stdn']:.3f}m σE={r['stde']:.3f}m "
                    f"σU={r['stdu']:.3f}m")
            pms.append(
                f'    <Placemark><description>{desc}</description>'
                f'<styleUrl>#{dot_style}</styleUrl>'
                f'<Point><altitudeMode>absolute</altitudeMode>'
                f'<coordinates>{r["lon"]},{r["lat"]},{r["hgt"]}</coordinates>'
                f'</Point></Placemark>')
        return '\n'.join(pms)

    kml = f'''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>PocketSDR Position{tag}</name>
  <description>Altitude reference: {alt_note}</description>
  <Style id="pppLine"><LineStyle><color>ff0000ff</color><width>2</width></LineStyle></Style>
  <Style id="pppDot"><IconStyle><color>ff0000ff</color><scale>0.5</scale></IconStyle>
                    <LabelStyle><scale>0</scale></LabelStyle></Style>
  <Style id="sppLine"><LineStyle><color>ff00ff00</color><width>1</width></LineStyle></Style>
  <Style id="sppDot"><IconStyle><color>ff00ff00</color><scale>0.4</scale></IconStyle>
                    <LabelStyle><scale>0</scale></LabelStyle></Style>
  <Folder>
    <name>PPP (Q=6) — {len(ppp)} epochs</name>
{track(ppp, "PPP", "pppLine", "pppDot")}
  </Folder>
  <Folder>
    <name>SPP (Q=5) — {len(spp)} epochs</name>
{track(spp, "SPP", "sppLine", "sppDot")}
  </Folder>
</Document>
</kml>'''
    with open(path, 'w') as f:
        f.write(kml)
    kind = f'KML{tag}:'
    print(f'{kind:<16} {path}  ({len(ppp)} PPP + {len(spp)} SPP)')

# ── GPX ───────────────────────────────────────────────────────────────────────

def write_gpx(rows, path, label=''):
    def trkpt(r):
        dt = gpst_to_utc(r)
        ts = dt.strftime('%Y-%m-%dT%H:%M:%S.') + f'{dt.microsecond//1000:03d}Z'
        ql = SOLQ_LABEL.get(r['q'], str(r['q']))
        return (f'      <trkpt lat="{r["lat"]:.9f}" lon="{r["lon"]:.9f}">'
                f'<ele>{r["hgt"]:.3f}</ele><time>{ts}</time>'
                f'<desc>{ql} ns={r["ns"]} σN={r["stdn"]:.3f}m</desc>'
                f'</trkpt>')

    ppp_pts = '\n'.join(trkpt(r) for r in rows if r['q'] == 6)
    spp_pts = '\n'.join(trkpt(r) for r in rows if r['q'] != 6)
    ppp_n = sum(1 for r in rows if r['q'] == 6)
    spp_n = len(rows) - ppp_n
    tag = f' ({label})' if label else ''

    gpx = f'''<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="PocketSDR pocket_export.py"
     xmlns="http://www.topografix.com/GPX/1/1">
  <trk>
    <name>PPP (Q=6){tag}</name>
    <trkseg>
{ppp_pts}
    </trkseg>
  </trk>
  <trk>
    <name>SPP (Q=5){tag}</name>
    <trkseg>
{spp_pts}
    </trkseg>
  </trk>
</gpx>'''
    with open(path, 'w') as f:
        f.write(gpx)
    kind = f'GPX{tag}:'
    print(f'{kind:<16} {path}  ({ppp_n} PPP + {spp_n} SPP)')

# ── GeoJSON ───────────────────────────────────────────────────────────────────

def write_geojson(rows, path, label=''):
    source = label or 'main'
    features = [{
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': [round(r['lon'], 9), round(r['lat'], 9),
                            round(r['hgt'], 3)],
        },
        'properties': {
            'source':         source,
            't_s':            r['t'],
            'datetime_gpst':  gpst_str(r),
            'quality':        r['q'],
            'quality_label':  SOLQ_LABEL.get(r['q'], str(r['q'])),
            'n_sats':         r['ns'],
            'stdn_m':         round(r['stdn'], 4),
            'stde_m':         round(r['stde'], 4),
            'stdu_m':         round(r['stdu'], 4),
            'dtr_s':          round(r['dtr'], 9),
        },
    } for r in rows]
    fc = {'type': 'FeatureCollection', 'features': features}
    with open(path, 'w') as f:
        json.dump(fc, f, separators=(',', ':'))
    tag = f' ({label})' if label else ''
    kind = f'GeoJSON{tag}:'
    print(f'{kind:<16} {path}  ({len(features)} features)')

# ── CSV / Excel ───────────────────────────────────────────────────────────────

def make_df(rows):
    """Build a DataFrame from $POS or $REFPOS row dicts.

    Optional fields (ecef, velocity, DOP, SPP seed clocks) use .get() so the
    same function works for both row types; missing values become NaN.
    SPP seed clock fields are merged into row dicts by main() before this call:
      spp_seed_stat / spp_seed_dtr_s  – main SPP seed (AOWR-corrected)
      ref_spp_stat  / ref_spp_dtr_s   – reference SPP seed (uncorrected, free clock)
    """
    _nan = float('nan')
    return pd.DataFrame([{
        't_s':           r['t'],
        'datetime_gpst': gpst_str(r),
        'gps_week':      (_wt := gpst_week_tow(r))[0],
        'tow_s':         _wt[1],
        # Geodetic
        'lat_deg':       r['lat'],
        'lon_deg':       r['lon'],
        'hgt_m':         r['hgt'],
        # ECEF position (NaN for $REFPOS rows which don't carry ECEF)
        'ecef_x_m':      r.get('ecef_x', _nan),
        'ecef_y_m':      r.get('ecef_y', _nan),
        'ecef_z_m':      r.get('ecef_z', _nan),
        # ECEF velocity (non-zero only in kinematic PPP; NaN for $REFPOS)
        'vel_x_m_s':     r.get('vel_x', _nan),
        'vel_y_m_s':     r.get('vel_y', _nan),
        'vel_z_m_s':     r.get('vel_z', _nan),
        # Solution quality
        'quality':       r['q'],
        'quality_label': SOLQ_LABEL.get(r['q'], str(r['q'])),
        'n_sats':        r['ns'],
        'stdn_m':        r['stdn'],
        'stde_m':        r['stde'],
        'stdu_m':        r['stdu'],
        'dtr_s':         r['dtr'],
        # Clock drift rate (s/s): Doppler-derived in SPP; WLS in SC AOWR; 0 in PPP
        'dtr_drift_s':   r.get('dtr_drift', _nan),
        # DOP (NaN for old logs without DOP fields)
        # For $POS (clock-corrected): gdop=0, pdop=PDOP3 (position-only, no clock coupling)
        # For $REFPOS (4-unknown):    gdop=GDOP, pdop=PDOP with clock-position coupling
        'gdop':          r.get('gdop', _nan),
        'pdop':          r.get('pdop', _nan),
        'hdop':          r.get('hdop', _nan),
        'vdop':          r.get('vdop', _nan),
        # SPP seed receiver clock offsets (debug; NaN when not in -ps_sc mode or no fix)
        'spp_seed_stat':  r.get('spp_seed_stat', _nan),
        'spp_seed_dtr_s': r.get('spp_seed_dtr',  _nan),
        'ref_spp_stat':   r.get('ref_spp_stat',  _nan),
        'ref_spp_dtr_s':  r.get('ref_spp_dtr',   _nan),
        # SC_AOWR clock estimate and drift rate (NaN when AOWR not yet active)
        'sc_aowr_clk_s':  r.get('sc_aowr_clk',   _nan),
        'sc_aowr_drift':  r.get('sc_aowr_drift',  _nan),
        # PPP KF convergence diagnostics from PPPOS_DBG (NaN when not in PPP mode
        # or old logs without these fields)
        'ppp_kf_pos_std_m': r.get('ppp_pos_std', _nan),
        'ppp_kf_ndual':     r.get('ppp_ndual',   _nan),
        'ppp_clk_m':        r.get('ppp_clk',     _nan),
        'ppp_clk_std_m':    r.get('ppp_clk_std', _nan),
        'ppp_zwd_m':        r.get('ppp_zwd',     _nan),
        'ppp_zwd_std_m':    r.get('ppp_zwd_std', _nan),
        'ppp_amb_std_m':    r.get('ppp_amb_std', _nan),
        'ppp_nconv':        r.get('ppp_nconv',   _nan),
        'ppp_res_c_m':      r.get('ppp_res_c',   _nan),
    } for r in rows])

def write_csv(rows, path, label=''):
    df = make_df(rows)
    df.to_csv(path, index=False, float_format='%.9g')
    tag = f' ({label})' if label else ''
    kind = f'CSV{tag}:'
    print(f'{kind:<16} {path}  ({len(df)} rows)')

def write_excel(rows, path, ref_rows=None):
    df = make_df(rows)
    if 'quality' in df.columns:
        ppp_df = df[df['quality'] == 6]
        spp_df = df[df['quality'] != 6]
    else:
        ppp_df = spp_df = df
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        df.to_excel(xw, sheet_name='All', index=False)
        ppp_df.to_excel(xw, sheet_name='PPP', index=False)
        spp_df.to_excel(xw, sheet_name='SPP', index=False)
        if ref_rows:
            rdf = make_df(ref_rows)
            rppp = rdf[rdf['quality'] == 6]
            rspp = rdf[rdf['quality'] != 6]
            rdf.to_excel(xw, sheet_name='Ref', index=False)
            rppp.to_excel(xw, sheet_name='RefPPP', index=False)
            rspp.to_excel(xw, sheet_name='RefSPP', index=False)
    ref_note = f' + {len(ref_rows)} ref' if ref_rows else ''
    print(f'Excel:           {path}  ({len(df)} rows{ref_note})')

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Export PocketSDR log to KML / GPX / GeoJSON / CSV / Excel.\n'
                    '$REFPOS records (SC AOWR reference solver) are exported '
                    'automatically with --all, or via --ref-* flags.')
    ap.add_argument('logfile', help='pocket.log file')
    ap.add_argument('--all', action='store_true',
                    help='write all formats; names derived from logfile stem')
    ap.add_argument('--ellipsoidal-height', action='store_true',
                    help='KML altitude = raw WGS-84 ellipsoidal height '
                         '(default: convert to EGM96/MSL height via pyproj, '
                         "matching Google Earth's absolute altitude mode)")
    # Main ($POS) output
    ap.add_argument('--kml',     metavar='FILE')
    ap.add_argument('--gpx',     metavar='FILE')
    ap.add_argument('--geojson', metavar='FILE')
    ap.add_argument('--csv',     metavar='FILE')
    ap.add_argument('--excel',   metavar='FILE')
    # Reference ($REFPOS) output
    ap.add_argument('--ref-kml',     metavar='FILE')
    ap.add_argument('--ref-gpx',     metavar='FILE')
    ap.add_argument('--ref-geojson', metavar='FILE')
    ap.add_argument('--ref-csv',     metavar='FILE')
    ap.add_argument('--ref-pos-file', metavar='FILE',
                    help='use an offline rnx2rtkp .pos file (e.g. standalone '
                         'PPP post-processing) as the reference-solver data '
                         'instead of pocket.log\'s own $REFPOS records -- for '
                         'sessions where the onboard reference solver got '
                         'stuck (see AOWR PPP KF reset bug in sdr_pvt.c)')
    ap.add_argument('--ref-dtr-debug', metavar='FILE',
                    help='DTRDBG debug trace (stderr from a temporarily-'
                         'instrumented rnx2rtkp -- see parse_dtr_debug_log() '
                         'docstring) giving true per-epoch receiver clock '
                         'bias for --ref-pos-file rows. Without this, dtr is '
                         'only approximated from the .pos timestamp to ~1 ms '
                         'precision -- use this flag when accurate dtr_s '
                         'matters (e.g. cross-checking against the main '
                         "solver's own clock estimate).")
    args = ap.parse_args()

    rows = parse_log(args.logfile)
    if rows:
        ppp_n = sum(1 for r in rows if r['q'] == 6)
        spp_n = len(rows) - ppp_n
        print(f'Parsed {len(rows)} $POS records from {args.logfile}'
              f'  ({ppp_n} PPP, {spp_n} SPP)')
    else:
        print(f'No $POS records in {args.logfile}', file=sys.stderr)

    if args.ref_pos_file:
        ref_rows = parse_rtklib_pos_as_ref(args.ref_pos_file, args.ref_dtr_debug)
        if ref_rows:
            rppp_n = sum(1 for r in ref_rows if r['q'] == 6)
            rspp_n = len(ref_rows) - rppp_n
            print(f'Parsed {len(ref_rows)} rows from {args.ref_pos_file} as '
                  f'reference solver data  ({rppp_n} PPP, {rspp_n} SPP) '
                  '-- overriding pocket.log $REFPOS')
    else:
        ref_rows = parse_refpos_log(args.logfile)
        if ref_rows:
            rppp_n = sum(1 for r in ref_rows if r['q'] == 6)
            rspp_n = len(ref_rows) - rppp_n
            print(f'Parsed {len(ref_rows)} $REFPOS records'
                  f'  ({rppp_n} PPP, {rspp_n} SPP)')

    if not rows and not ref_rows:
        print('No position records found, exiting.', file=sys.stderr)
        sys.exit(1)

    # Merge SPP seed clock entries into rows for debug export.
    # REF_SPP_SEED (PPP mode) and REF_SPP (SPP mode) both carry the reference
    # solver clock; combine them into a single lookup.
    spp_clk = parse_spp_seed_log(args.logfile, 'SPP_SEED')
    ref_spp_clk = {**parse_spp_seed_log(args.logfile, 'REF_SPP_SEED'),
                   **parse_spp_seed_log(args.logfile, 'REF_SPP')}
    if spp_clk or ref_spp_clk:
        for r in rows:
            if r['t'] in spp_clk:
                r['spp_seed_stat'] = spp_clk[r['t']]['stat']
                r['spp_seed_dtr']  = spp_clk[r['t']]['dtr']
            if r['t'] in ref_spp_clk:
                r['ref_spp_stat'] = ref_spp_clk[r['t']]['stat']
                r['ref_spp_dtr']  = ref_spp_clk[r['t']]['dtr']
        print(f'Merged {len(spp_clk)} SPP_SEED and {len(ref_spp_clk)} REF_SPP clock entries')

    sc_aowr = parse_sc_aowr_log(args.logfile)
    if sc_aowr:
        n_reconstructed = 0
        for r in rows:
            if r['t'] in sc_aowr:
                r['sc_aowr_clk']   = sc_aowr[r['t']]['clk']
                r['sc_aowr_drift'] = sc_aowr[r['t']]['drift']
                # Logs from before the out_log_pos() CLIGHT-division fix: in
                # AOWR clock-fixed mode, dtr[0] is already seconds (retained
                # from clock_est), but out_log_pos() divided it by CLIGHT
                # anyway whenever stat==SOLQ_PPP, landing at ~1e-11 s which
                # rounds to exactly 0.000000000 at the log's %.9f precision.
                # sc_aowr_clk is that same clock_est value (same 'time' tick,
                # written by the same update_sol() call) -- use it to recover
                # dtr_s for exactly the rows showing that corruption
                # signature (q=6 PPP, dtr logged as exactly zero).
                if r['q'] == 6 and r['dtr'] == 0.0:
                    r['dtr'] = sc_aowr[r['t']]['clk']
                    n_reconstructed += 1
        print(f'Merged {len(sc_aowr)} SC_AOWR clock/drift entries')
        if n_reconstructed:
            print(f'Reconstructed dtr_s for {n_reconstructed} PPP epoch(s) '
                  '(pre-fix CLIGHT-division bug in out_log_pos(); see '
                  'sc_aowr_clk_s for the source value)')

    pppos_dbg = parse_pppos_dbg_log(args.logfile)
    if pppos_dbg:
        for r in rows:
            if r['t'] in pppos_dbg:
                d = pppos_dbg[r['t']]
                r['ppp_pos_std'] = d['pos_std']
                r['ppp_ndual']   = d['ndual']
                r['ppp_clk']     = d['clk']
                r['ppp_clk_std'] = d['clk_std']
                r['ppp_zwd']     = d['zwd']
                r['ppp_zwd_std'] = d['zwd_std']
                r['ppp_amb_std'] = d['amb_std']
                r['ppp_nconv']   = d['nconv']
                r['ppp_res_c']   = d['res_c']
        print(f'Merged {len(pppos_dbg)} PPPOS_DBG convergence entries')

    stem = os.path.splitext(args.logfile)[0]
    ref_stem = stem + '_ref'

    want = args.all or any([args.kml, args.gpx, args.geojson, args.csv, args.excel,
                            args.ref_kml, args.ref_gpx, args.ref_geojson, args.ref_csv])
    if not want:
        ap.print_help()
        sys.exit(0)

    use_geoid = not args.ellipsoidal_height

    # ── main ($POS) output ──
    if args.all or args.kml:
        write_kml(rows,     args.kml     or stem + '.kml', use_geoid=use_geoid)
    if args.all or args.gpx:
        write_gpx(rows,     args.gpx     or stem + '.gpx')
    if args.all or args.geojson:
        write_geojson(rows, args.geojson or stem + '.geojson')
    if args.all or args.csv:
        write_csv(rows,     args.csv     or stem + '.csv')
    if args.all or args.excel:
        # Include ref sheets in the same Excel file for easy comparison
        write_excel(rows,   args.excel   or stem + '.xlsx',
                    ref_rows=ref_rows or None)

    # ── reference ($REFPOS) output ──
    if ref_rows:
        if args.all or args.ref_kml:
            write_kml(ref_rows,     args.ref_kml     or ref_stem + '.kml',
                      label='ref', use_geoid=use_geoid)
        if args.all or args.ref_gpx:
            write_gpx(ref_rows,     args.ref_gpx     or ref_stem + '.gpx',
                      label='ref')
        if args.all or args.ref_geojson:
            write_geojson(ref_rows, args.ref_geojson or ref_stem + '.geojson',
                          label='ref')
        if args.all or args.ref_csv:
            write_csv(ref_rows,     args.ref_csv     or ref_stem + '.csv',
                      label='ref')
    elif any([args.ref_kml, args.ref_gpx, args.ref_geojson, args.ref_csv]):
        print('Warning: no $REFPOS records found — ref output skipped',
              file=sys.stderr)

if __name__ == '__main__':
    main()
