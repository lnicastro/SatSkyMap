#!/bin/bash
# 
# Can also retrive TLE data from Space-Track (https://www.space-track.org/).
#
# This script uses CelesTrak (https://celestrak.com/NORAD/elements/)
#
# To get TLE for a single sat, use the NORAD Id, e.g.: https://celestrak.com/satcat/tle.php?CATNR=44874
#
# Edit to match your needs.
#
# LN @ INAF-OAS Jan. 2020.  Last change: 05/05/2021
#--

set +o noclobber

URL="https://celestrak.com/NORAD/elements/"

TIME="10"

LOGFILE=/tmp/tle_retrieve.log
BINDIR=/usr/local/bin
TMPFILE=tle.tmp

# Where to store the files
if [ $# -eq 0 ]; then
	OUTDIR=/usr/local/TLErepo
else
	OUTDIR=$1
fi

# Goto destination directory
cd $OUTDIR

# List of most relevant TLE files
tles=( tle-new stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos planet spire geo intelsat ses iridium iridium-NEXT starlink orbcomm globalstar amateur x-comm other-comm satnogs gorizont raduga molniya  gps-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other )

date | tee $LOGFILE

for TLE_NAME in "${tles[@]}"
do
	echo -n "===> $TLE_NAME ..." | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	#curl -s -O --max-time $TIME "$URL$TLE"
	if [ `curl -s -w '%{http_code}' -o $TMPFILE --max-time $TIME "$URL$TLE"` -eq 200 ]; then  # Check for success code
		mv $TMPFILE $TLE
		echo retrieved. | tee -a $LOGFILE
	else
		echo retrieval failed. | tee -a $LOGFILE
	fi

done


# List of Debris TLE files
tlesdeb=( 2019-006 1999-025 iridium-33-debris cosmos-2251-debris )

for TLE_NAME in "${tlesdeb[@]}"
do
	echo -n "===> $TLE_NAME ..." | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	if [ `curl -s -w '%{http_code}' -o $TMPFILE --max-time $TIME "$URL$TLE"` -eq 200 ]; then
		mv $TMPFILE $TLE
		echo retrieved. | tee -a $LOGFILE
	else
		echo Retrieval failed. | tee -a $LOGFILE
	fi

done


# TLE files with different url/extension name (in Communications Satellites section)
  echo -n "===> /satcat/gpz.php ..." | tee -a $LOGFILE
  if [ `curl -s -w '%{http_code}' -o $TMPFILE --max-time $TIME https://celestrak.com/satcat/gpz.php` -eq 200 ]; then
	mv $TMPFILE gpz.txt
	echo retrieved. | tee -a $LOGFILE
  else
	echo Retrieval failed. | tee -a $LOGFILE
  fi

  echo -n "===> /satcat/gpz-plus.php ..." | tee -a $LOGFILE
  if [ `curl -s -w '%{http_code}' -o $TMPFILE --max-time $TIME https://celestrak.com/satcat/gpz-plus.php` -eq 200 ]; then
	mv $TMPFILE gpz-plus.txt
	echo retrieved. | tee -a $LOGFILE
  else
	echo Retrieval failed. | tee -a $LOGFILE
  fi


#
# Produce the TLE list file for special satellites. Comment out if not needed.
#
awk '/HST|CXO|GLAST|SWIFT|NUSTAR|AGILE|INTEGRAL|ASTROSAT|HXMT|XMM|WISE|SDO/ { print ; for(n=0; n<2; n++) { getline ; print } }' science.txt | tr -d '\r' > special.txt
awk '/ZARYA/ { print ; for(n=0; n<2; n++) { getline ; print } }' stations.txt | tr -d '\r' >> special.txt
awk '/CHEOPS|TESS/ { print ; for(n=0; n<2; n++) { getline ; print } }' active.txt | tr -d '\r' >> special.txt


# Add Gaia TLE from 2021 list (see https://github.com/Bill-Gray/tles/ => 13074a21.tle)
  if [ ! -f $OUTDIR/gaia_2021.txt ]; then
	echo To have Gaia TLEs move gaia_2021.txt from the scripts to the TLE dir.
  else
#MJD=`$BINDIR/mjdnow.php | sed -e 's/[^\.]*$//'`  # Remove fractional part
# Use approx int MJD. See also the C code in src directory.
	MJD=`expr $(date +%s) / 86400 + 40588`
	awk -v mjd="${MJD}" '$0 ~ mjd { print ; for(n=0; n<3; n++) { getline ; if ( index($1, "Gaia") > 0 ) { print "Gaia"} else { print } } }' gaia_2021.txt >> special.txt
  fi


#
# Produce the updated list of number of sats in the TLE file. Comment out if not needed.
#
$BINDIR/tle_satcount.sh

date | tee -a $LOGFILE
