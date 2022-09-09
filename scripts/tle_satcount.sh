#!/bin/bash
# 
# Count and log the number of satellites from the SpaceTrack retrieved TLE data
#
# Edit to match your needs.
#
#
# LN @ INAF-OAS Jan. 2020.  Last change: 09/09/2022
#--

set +o noclobber

#NOW=$(date +"%Y-%m-%d-%H%M")
#NOW=$(date +"%Y%m%d")

# Output log file. First tw in /tmp the list of all current NORAD IDs in TLE dir.
#/tmp/tle_satcount_$NOW.log
OUTFILE=/tmp/tle_satcount.log
OUTFILE_DEB=/tmp/tle_satcount_deb.log
OUTFILE_NORAD=cur_norad_ids.lis

# Where the TLE files are stored
if [ $# -eq 0 ]; then
	TLEDIR=/usr/local/TLErepo
else
	TLEDIR=$1
fi

# Goto TLE dir
cd $TLEDIR

# Remove old log files (if they exist)
rm -f $OUTFILE $OUTFILE_DEB $OUTFILE_NORAD

# List of most relevant TLE files, excluding debris
tlesmain=( last-30-days stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos spire geo iridium iridium-NEXT globalstar swarm amateur x-comm other-comm satnogs gorizont raduga molniya gnss gps-ops glo-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other gpz gpz-plus )

# Supplemental TLE files not in main list
# Note: "iss" not includede and "gps" is renamed to "gps-ops" => removed from the list
tlessupp=( starlink oneweb planet glonass meteosat intelsat ses telesat orbcomm cpf )

# List of debris TLE files
tlesdeb=( 1982-092 1999-025 iridium-33-debris cosmos-2251-debris )

for TLE_NAME in "${tlesmain[@]}"
do
	TLE=$TLE_NAME".txt"
	wc -l $TLE | awk '{print substr($2, 0, index($2, ".")-1)" "$1/3}' | tee -a $OUTFILE

	sed -n 3~3p $TLE | cut -d " " -f 2 >> ${OUTFILE_NORAD}_tmp
done

for TLE_NAME in "${tlessupp[@]}"
do
	TLE=$TLE_NAME".txt"
	wc -l $TLE | awk '{print substr($2, 0, index($2, ".")-1)" "$1/3}' | tee -a $OUTFILE

	sed -n 3~3p $TLE | cut -d " " -f 2 >> ${OUTFILE_NORAD}_tmp
done

# Debris names are all the same in a file
for TLE_NAME in "${tlesdeb[@]}"
do
	TLE=$TLE_NAME".txt"
	wc -l $TLE | awk '{print substr($2, 0, index($2, ".")-1)" "$1/3}' | tee -a $OUTFILE_DEB
done

# Count "unique" NORAD IDs
sort -u < ${OUTFILE_NORAD}_tmp > $OUTFILE_NORAD
rm ${OUTFILE_NORAD}_tmp

n_sats=`cat $OUTFILE_NORAD | wc -l`

# Debris IDs are all different so just sum the counts for each file
n_debs=`awk '{sum += $2} END {print sum}' $OUTFILE_DEB`

# The total
n_tot=$((n_sats + n_debs))

echo ""
echo "$OUTFILE | $OUTFILE_DEB | $OUTFILE_NORAD ready"

#
# Create a PHP array. Comment if not needed.
#
PHPOUT=cur_tlesat_count.php

echo '<?php
$my_satcount_arr = array(' > $PHPOUT

awk '{print "\""$1"\"" " => "$2","}' $OUTFILE >> $PHPOUT
awk '{print "\""$1"\"" " => "$2","}' $OUTFILE_DEB >> $PHPOUT

# Here we add the count of (unique) satellites including / excluding debris
echo "'@ALL_merged_nodeb' => $n_sats," >> $PHPOUT
echo "'@ALL_merged' => $n_tot" >> $PHPOUT

echo ');
?>' >> $PHPOUT

echo $PHPOUT ready
