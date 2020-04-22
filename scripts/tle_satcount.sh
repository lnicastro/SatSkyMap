#!/bin/bash
# 
# Count and log the number of satellites from the SpaceTrack retrived TLE data
#
# Edit to match your needs.
#
# LN @ INAF-OAS Jan. 2020.  Last change: 22/04/2020
#--

set +o noclobber

#NOW=$(date +"%Y-%m-%d-%H%M")
NOW=$(date +"%Y%m%d")

# Output log file
OUTFILE=/tmp/tle_satcount_$NOW.log

# Where the TLE files are stored
if [ $# -eq 0 ]; then
	TLEDIR=/usr/local/TLErepo
else
	TLEDIR=$1
fi

# Goto TLE dir
cd $TLEDIR

# Remove old log file (if it exists)
rm -f $OUTFILE


# List of most relevant TLE files, including debris
tles=( tle-new stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos planet spire geo intelsat ses iridium iridium-NEXT starlink orbcomm globalstar amateur x-comm other-comm satnogs gorizont raduga molniya  gps-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other 2019-006 1999-025 iridium-33-debris cosmos-2251-debris gpz.txt gpz-plus.txt )

for TLE_NAME in "${tles[@]}"
do
	TLE=$TLE_NAME".txt"
	wc -l $TLE | awk '{print substr($2, 0, index($2, ".")-1)" "$1/3}' | tee -a $OUTFILE

done

echo ""
echo "$OUTFILE ready"

#
# Create a PHP array. Comment if not needed.
#
PHPOUT=cur_tlesat_count.php

echo '<?php
$my_satcount_arr = array(' > $PHPOUT

awk '{print "\""$1"\"" " => "$2","}' $OUTFILE >> $PHPOUT

# Should remove last ,
echo ');
?>' >> $PHPOUT

echo $PHPOUT ready
