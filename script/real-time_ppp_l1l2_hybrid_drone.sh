#! /bin/bash
#

trap 'stty sane' INT TERM

###
# get absolute path to conf dir and sibling repos
scriptdir=$(cd $(dirname $0) && pwd)
confpath=$scriptdir/../conf
pythondir=$scriptdir/../python
conffile=$confpath/pocket_L1L2_4MHz_ext10MHz_ch34.conf
aowr_file=$scriptdir/../../jaxa-asyncOWR-prototype/work/dt_aowr_gnss.txt

# Optional: precise ephemeris files for faster PPP convergence.
# Download IGS ultra-rapid products before the experiment:
#   wget https://igs.ign.fr/pub/igs/products/ultra/WEEK/DDDH/IGU<WEEK><DOW>_<HH>.SP3.Z
#   wget https://igs.ign.fr/pub/igs/products/ultra/WEEK/DDDH/IGU<WEEK><DOW>_<HH>.CLK_30S.Z
# Set SP3_FILE and CLK_FILE to the (uncompressed) paths, or leave empty.
SP3_FILE=""
CLK_FILE=""

build_nav_args() {
    local args=""
    [ -n "$SP3_FILE" ] && [ -f "$SP3_FILE" ] && args="$args -nav $SP3_FILE"
    [ -n "$CLK_FILE" ] && [ -f "$CLK_FILE" ] && args="$args -nav $CLK_FILE"
    echo "$args"
}

# Log folder
date_str="`date +'%Y%m%d%H%M%S'`"
abs_date_str="$(pwd)/$date_str"
mkdir "$abs_date_str"

cp $conffile $abs_date_str/
cp $confpath/ppp_kinematic_drone.ini $abs_date_str/
cp $0 $abs_date_str/
[ -n "$SP3_FILE" ] && cp "$SP3_FILE" $abs_date_str/ 2>/dev/null
[ -n "$CLK_FILE" ] && cp "$CLK_FILE" $abs_date_str/ 2>/dev/null

# Build nav args now (before writing runner script)
nav_args="$(build_nav_args)"

###
# Write runner script for the trk tmux window.
# Uses absolute paths so tmux can execute it from any working directory.
cat > "$abs_date_str/_run_trk.sh" << EOF
#! /bin/bash
trap 'stty sane' INT TERM

pocket_trk \\
    -c $conffile \\
    -sig L1CA,L2CM -prn 1-32,159,193-202 -ps_prn 159 -rfch 3,4 -f 4 -IQ 2 -ppp -ps_sc \\
    -ps_sc_file $aowr_file \\
    -opt $confpath/ppp_kinematic_drone.ini \\
    $nav_args \\
    -log $abs_date_str/pocket.log -nmea /dev/ttyUSB0 \\
    -debug $abs_date_str/trace.log

echo "[trk] pocket_trk exited -- running post-processing..."
cd "$abs_date_str"

python3 $pythondir/pocket_pos_plot.py pocket.log --out pos_result.pdf
python3 $pythondir/pocket_export.py pocket.log --all
python3 $pythondir/pocket_obs2rnx.py pocket.log

grep '^\$AOWR' pocket.log | \\
    awk -F',' 'BEGIN{print "time,year,month,day,hour,min,sec,prn,P,L,dt_raw,dt_pr,dt_cp,dt0,cp_thresh,count,outlier"} \\
    {print \$2","\$3","\$4","\$5","\$6","\$7","\$8","\$9","\$10","\$11","\$12","\$13","\$14","\$15","\$16","\$17","\$18}' \\
    > aowr.csv

echo "[trk] done. Log: $abs_date_str"
EOF
chmod +x "$abs_date_str/_run_trk.sh"

###
# Launch tmux session with two panes side by side
SESSION="pocketsdr"
tmux kill-session -t $SESSION 2>/dev/null
tmux new-session  -d -s $SESSION -n "drone"

# Split into left (trk) and right (aowr) panes
tmux split-window -h -t "$SESSION:drone"

# Left pane (0): run pocket_trk then post-processing automatically
tmux send-keys -t "$SESSION:drone.0" "bash $abs_date_str/_run_trk.sh" Enter

# Right pane (1): cd to the aowr work dir, ready for the user to start ranging
tmux send-keys -t "$SESSION:drone.1" "cd $scriptdir/../../jaxa-asyncOWR-prototype/work" Enter

# Focus on left pane (trk)
tmux select-pane -t "$SESSION:drone.0"
tmux attach-session -t $SESSION
