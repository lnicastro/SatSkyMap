#!/bin/bash
# 
# Retrive TLE data from SpaceTrack
#
# See: https://celestrak.com/NORAD/elements/
#
# LN @ INAF-OAS Jan. 2020.  Last change: 18/02/2020
#--

set +o noclobber

URL="https://celestrak.com/NORAD/elements/"

TIME="10"

LOGFILE=/tmp/tle_retrieve.log

# Where to store the files
if [ $# -eq 0 ]; then
	OUTDIR=/usr/local/TLErepo
else
	OUTDIR=$1
fi

cd $OUTDIR

# List of most relevant TLE files
tles=( tle-new stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos planet spire geo intelsat ses iridium iridium-NEXT starlink orbcomm globalstar amateur x-comm other-comm satnogs gorizont raduga molniya  gps-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other )

date | tee $LOGFILE

for TLE_NAME in "${tles[@]}"
do
	echo "===> $TLE_NAME" | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	curl -s -O --max-time $TIME "$URL$TLE"

done


# List of Debries TLE files
tles=( 2019-006 1999-025 iridium-33-debris cosmos-2251-debris )

for TLE_NAME in "${tles[@]}"
do
	echo "===> $TLE_NAME" | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	curl -s -O --max-time $TIME "$URL$TLE"

done


# TLE files with different url/extension name (in Communications Satellites section)
	echo "===> /satcat/gpz.php" | tee -a $LOGFILE
	curl -s -o gpz.txt https://celestrak.com/satcat/gpz.php

	echo "===> /stacat/gpz-plus.php" | tee -a $LOGFILE
	curl -s -o gpz-plus.txt https://celestrak.com/satcat/gpz-plus.php

date | tee -a $LOGFILE
