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

trap 'stty sane' INT TERM

###
# get absolute path to conf dir
confpath=$(cd $(dirname $0) && pwd)/../conf
conffile=$confpath/pocket_L1L2_4MHz_ext10MHz_ch34.conf

# Log folder
date_str="`date +'%Y%m%d%H%M%S'`"
mkdir $date_str

cp $conffile $date_str

###
# run
pocket_trk \
    -c $conffile \
    -sig L1CA -prn 1-32,193-202 -rfch 3 -f 4 -IQ 2 -pppsg \
    -log $date_str/pocket.log -nmea /dev/ttyUSB0 \
    -debug $date_str/trace.log

###
# plot and save position figure
pushd .
cd $date_str
python3 ../../python/pocket_pos_plot.py pocket.log --out pos_result.pdf
popd
