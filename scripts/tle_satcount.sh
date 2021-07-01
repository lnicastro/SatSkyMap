#!/bin/bash
# 
# Count and log the number of satellites from the SpaceTrack retrieved TLE data
#
# Edit to match your needs. On MacoS "sed" needs to be raplaced with "gsed" (gnu-sed) to work.
#
# LN @ INAF-OAS Jan. 2020.  Last change: 01/07/2021
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
#tlesmain=( tle-new stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos planet spire geo intelsat ses iridium iridium-NEXT starlink orbcomm globalstar amateur x-comm other-comm satnogs gorizont raduga molniya gps-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other gpz gpz-plus )
tlesmain=( tle-new stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos spire geo iridium iridium-NEXT globalstar swarm amateur x-comm other-comm satnogs gorizont raduga molniya gnss gps-ops glo-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other gpz gpz-plus )

# Supplemental TLE files (not that gps is renamed to gps-ops)
tlessupp=( glonass intelsat oneweb orbcomm planet ses starlink transporter-2 meteosat telesat iss cpf )

# List of Debris TLE files
tlesdeb=( 2019-006 1999-025 iridium-33-debris cosmos-2251-debris )

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

# Debries names are all the same in a file
for TLE_NAME in "${tlesdeb[@]}"
do
	TLE=$TLE_NAME".txt"
	wc -l $TLE | awk '{print substr($2, 0, index($2, ".")-1)" "$1/3}' | tee -a $OUTFILE_DEB
done

# Count "unique" NORAD IDs
sort -u < ${OUTFILE_NORAD}_tmp > $OUTFILE_NORAD
rm ${OUTFILE_NORAD}_tmp

n_sats=`cat $OUTFILE_NORAD | wc -l`

# Debries IDs are all different so just sum the counts for each file
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

# Here we add the count of (unique) satellites including / excluding debris
echo "'@ALL_merged_nodeb' => $n_sats," >> $PHPOUT
echo "'@ALL_merged' => $n_tot" >> $PHPOUT

echo ');
?>' >> $PHPOUT

echo $PHPOUT ready
