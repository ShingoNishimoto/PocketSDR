#! /bin/bash
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
    -sig L1CA,L2CM -prn 1-32,193-202 -rfch 3,4 -f 4 -IQ 2 -ppps \
    -opt $confpath/ppp_l2_opt.ini \
    -log $date_str/pocket.log -nmea /dev/ttyUSB0 \
    -debug $date_str/trace.log

###
# plot and save position figure
pushd .
cd $date_str
python3 ../../python/pocket_pos_plot.py pocket.log --out pos_result.pdf
popd