#!/bin/bash
#
# post_process_rtk.sh - full offline RTK post-processing + comparison
# pipeline for one pocket_trk session:
#   1. rover.obs/rover.nav (offline) + base obs fetch + rnx2rtkp
#        [pocket_rtk_postproc.py]
#   2. rtk.pos trajectory sanity-check plot
#        [pocket_pos_plot.py]
#   3. onboard main/reference solver vs RTK-truth comparison
#        [pocket_rtk_compare.py]
#   4. rtk.pos full position/state history -> Excel + KML
#        [pocket_rtk_export.py]
#
# Usage:
#   script/post_process_rtk.sh <session_dir> [extra pocket_rtk_postproc.py args...]
#
# Examples:
#   script/post_process_rtk.sh work/20260705162229
#   script/post_process_rtk.sh work/20260705162229 --skip-fetch \
#       --rnx2rtkp ../RTKLIB_fork/app/consapp/rnx2rtkp/gcc/rnx2rtkp
#
trap 'stty sane' INT TERM

if [ -z "$1" ]; then
    echo "Usage: $0 <session_dir> [extra pocket_rtk_postproc.py args...]"
    exit 1
fi
session_dir=${1%/}
shift

scriptdir=$(cd "$(dirname "$0")" && pwd)
pydir=$scriptdir/../python

echo "=== 1/3: RTK post-processing ==="
python3 "$pydir/pocket_rtk_postproc.py" "$session_dir" "$@"
if [ ! -s "$session_dir/rtk.pos" ]; then
    echo "error: $session_dir/rtk.pos was not produced -- aborting" >&2
    exit 1
fi

echo
echo "=== 2/3: rtk.pos sanity-check plot ==="
python3 "$pydir/pocket_pos_plot.py" "$session_dir/rtk.pos" \
    --out "$session_dir/rtk_check.png" \
    || echo "warning: plotting rtk.pos failed (see above); continuing" >&2

echo
echo "=== 3/4: main/reference solver vs RTK truth ==="
python3 "$pydir/pocket_rtk_compare.py" "$session_dir/pocket.log" "$session_dir/rtk.pos" \
    --out "$session_dir/aowr_compare.png" --csv "$session_dir/aowr_compare.csv" \
    || echo "warning: pocket_rtk_compare.py failed (see above)" >&2

echo
echo "=== 4/4: rtk.pos position/state history -> Excel + KML ==="
python3 "$pydir/pocket_rtk_export.py" "$session_dir/rtk.pos" \
    --excel "$session_dir/rtk.xlsx" \
    --kml "$session_dir/rtk.kml" \
    || echo "warning: pocket_rtk_export.py failed (see above)" >&2

echo
echo "Done. Outputs in $session_dir/:"
echo "  rover.obs, rover.nav   - offline RINEX from pocket.log"
echo "  rtk.pos                - RTK post-processed trajectory"
echo "  rtk_check.png          - trajectory sanity-check plot"
echo "  aowr_compare.png/.csv  - main vs reference solver vs RTK truth"
echo "  rtk.xlsx               - full rtk.pos position/state history"
echo "  rtk.kml                - rtk.pos trajectory for Google Earth"
