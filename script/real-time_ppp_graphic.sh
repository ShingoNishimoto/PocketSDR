#! /bin/bash
#
# Single-frequency PPP using GRAPHIC (GRadient And PHase Ionospheric Correction)
# combination: G = (P_L1 + L_L1) / 2
#
# The GRAPHIC combination cancels first-order ionospheric delay using L1CA only.
# Tropospheric delay is estimated by the PPP Kalman filter (TROPOPT_ESTG).
# Requires: CH3 (L1CA, rfch 3) only — CH1 must run for its CLKOUT to drive CH3.
#
# Usage: ./real-time_ppp_graphic.sh
#
# Optional: set SP3_FILE and/or CLK_FILE to precise ephemeris paths for better
# orbit/clock accuracy. Download IGS ultra-rapid products before the experiment:
#   wget https://igs.ign.fr/pub/igs/products/ultra/WEEK/DDDH/IGU<WEEK><DOW>_<HH>.SP3.Z
#   wget https://igs.ign.fr/pub/igs/products/ultra/WEEK/DDDH/IGU<WEEK><DOW>_<HH>.CLK_30S.Z
#

trap 'stty sane' INT TERM

###
# get absolute path to conf dir
confpath=$(cd $(dirname $0) && pwd)/../conf
conffile=$confpath/pocket_L1L2_4MHz_ext10MHz_ch34.conf

# Optional: precise ephemeris files for better orbit/clock accuracy.
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
mkdir $date_str

cp $conffile $date_str
[ -n "$SP3_FILE" ] && cp "$SP3_FILE" $date_str/ 2>/dev/null
[ -n "$CLK_FILE" ] && cp "$CLK_FILE" $date_str/ 2>/dev/null

###
# run
pocket_trk \
    -c $conffile \
    -sig L1CA -prn 1-32,193-202 -rfch 3 -f 4 -IQ 2 -pppsg \
    -opt $confpath/ppp_graphic_opt.ini \
    $(build_nav_args) \
    -log $date_str/pocket.log -nmea /dev/ttyUSB0 \
    -debug $date_str/trace.log

###
# plot and export position data
pushd .
cd $date_str
python3 ../../python/pocket_pos_plot.py pocket.log --out pos_result.pdf
python3 ../../python/pocket_export.py pocket.log --all
python3 ../../python/pocket_obs2rnx.py pocket.log --sys G,J
popd
