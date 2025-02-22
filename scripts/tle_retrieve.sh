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
# 18/08/2022: Updated celestrak retrieval links
# 06/09/2022: Updated for iridium -> iridium-NEXT in the supplemental archive
# 13/09/2022: Added IXPE to the "special" TLE collection
# 05/07/2023: Added EUCLID to the "special" TLE collection
# 08/07/2023: Added JUICE to the "special" TLE collection
# 22/04/2024: Updated Gaia, Euclid, JWT, JUICE TLEs
# 15/01/2025: Updated Gaia, Euclid, JWT, JUICE (2024) TLEs
#
#
# LN @ INAF-OAS Jan. 2020.  Last change: 15/01/2025
#--

set +o noclobber

URL="https://celestrak.com/NORAD/elements/gp.php?FORMAT=tle&GROUP="
URLS="https://celestrak.com/NORAD/elements/gp.php?FORMAT=tle&SPECIAL="
URLSUPP="https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FORMAT=tle&FILE="

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
tlesmain=( last-30-days stations visual active analyst weather noaa goes resource sarsat dmc tdrss argos spire geo iridium iridium-NEXT globalstar swarm amateur x-comm other-comm satnogs gorizont raduga molniya gnss gps-ops glo-ops galileo beidou sbas nnss musson science geodetic engineering education military radar cubesat other )

date | tee $LOGFILE

for TLE_NAME in "${tlesmain[@]}"
do
	echo -n "===> $TLE_NAME ..." | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	if [ `curl -L -s -w '%{http_code}' -o $TMPFILE --max-time $TIME "$URL$TLE_NAME"` -eq 200 ]; then  # Check for success code
		mv $TMPFILE $TLE
		echo retrieved. | tee -a $LOGFILE
	else
		echo retrieval failed. | tee -a $LOGFILE
	fi

done

# Supplemental TLE files (see https://celestrak.com/NORAD/elements/supplemental/).
# Previewsly in the main list: gps (as "gps-ops"), intelsat, oneweb, orbcomm, planet, ses, starlink
# Note: multiple entries marking pre and post-maneuver [PM] TLEs (e.g. for intelsat) are not yet managed.
#       iss.txt not used (it reports "segments").
#
# Note: iridium is in fact iridium-NEXT

tlessupp=( starlink oneweb planet iridium gps glonass meteosat intelsat ses telesat orbcomm iss cpf )

echo ""
echo Supplemental GP Element Sets
echo ""

for TLE_NAME in "${tlessupp[@]}"
do
	echo -n "===> $TLE_NAME ..." | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	if [ `curl -L -s -w '%{http_code}' -o $TMPFILE --max-time $TIME $URLSUPP$TLE_NAME` -eq 200 ]; then
		if [ $TLE_NAME = iridium ]; then
			mv $TMPFILE iridium-NEXT.txt
		else
			mv $TMPFILE $TLE
		fi

		echo retrieved. | tee -a $LOGFILE
	else
		echo Retrieval failed. | tee -a $LOGFILE
	fi

done

# Patch for gps file name mismatch
if [ -f gps.txt ]; then
	mv gps.txt gps-ops.txt
fi



# List of Debris TLE files
tlesdeb=( 1982-092 1999-025 iridium-33-debris cosmos-2251-debris )  # 2019-006 replaced with 1982-092

for TLE_NAME in "${tlesdeb[@]}"
do
	echo -n "===> $TLE_NAME ..." | tee -a $LOGFILE

	TLE=$TLE_NAME".txt"
	if [ `curl -L -s -w '%{http_code}' -o $TMPFILE --max-time $TIME "$URL$TLE_NAME"` -eq 200 ]; then
		mv $TMPFILE $TLE
		echo retrieved. | tee -a $LOGFILE
	else
		echo Retrieval failed. | tee -a $LOGFILE
	fi

done


# TLE files with different url/extension name (in Communications Satellites section)
  echo -n "===> SPECIAL gpz ..." | tee -a $LOGFILE

  if [ `curl -L -s -w '%{http_code}' -o $TMPFILE --max-time $TIME "${URLS}gpz"` -eq 200 ]; then
	mv $TMPFILE gpz.txt
	echo retrieved. | tee -a $LOGFILE
  else
	echo Retrieval failed. | tee -a $LOGFILE
  fi

  echo -n "===> SPECIAL gpz-plus ..." | tee -a $LOGFILE
  if [ `curl -L -s -w '%{http_code}' -o $TMPFILE --max-time $TIME "${URLS}gpz-plus"` -eq 200 ]; then
	mv $TMPFILE gpz-plus.txt
	echo retrieved. | tee -a $LOGFILE
  else
	echo Retrieval failed. | tee -a $LOGFILE
  fi


#
# Produce the TLE list file for special satellites. Comment out if not needed.
#
awk '/HST|CXO|GLAST|SWIFT|NUSTAR|AGILE|INTEGRAL|ASTROSAT|HXMT|XMM|WISE|SDO|IXPE/ { print ; for(n=0; n<2; n++) { getline ; print } }' science.txt | tr -d '\r' > special.txt
awk '/ZARYA/ { print ; for(n=0; n<2; n++) { getline ; print } }' stations.txt | tr -d '\r' >> special.txt
awk '/CHEOPS|TESS/ { print ; for(n=0; n<2; n++) { getline ; print } }' active.txt | tr -d '\r' >> special.txt


# Add Gaia TLE from 2025 list (see https://github.com/Bill-Gray/tles/ => 13074a25.tle)
  if [ ! -f $OUTDIR/gaia_2025.txt ]; then
	echo To have Gaia TLEs move gaia_2025.txt from the scripts to the TLE dir.
  else
# Use approx int MJD. See also the C code in src directory.
	MJD=`expr $(date +%s) / 86400 + 40588`
	awk -v mjd="${MJD}" '$0 ~ mjd { print ; for(n=0; n<3; n++) { getline ; if ( index($1, "Gaia") > 0 ) { print "GAIA"} else { print } } }' gaia_2025.txt >> special.txt
  fi

# Add James Webb TLE from 2025 list (see https://github.com/Bill-Gray/tles/ => 21130a25.tle)
  if [ ! -f $OUTDIR/jwt_2025.txt ]; then
	echo To have JWT TLEs move jwt_2025.txt from the scripts to the TLE dir.
  else
# Use approx int MJD. See also the C code in src directory.
	MJD=`expr $(date +%s) / 86400 + 40588`
	awk -v mjd="${MJD}" '$0 ~ mjd { print ; for(n=0; n<3; n++) { getline ; if ( index($1, "James") > 0 ) { print "JAMES WEBB"} else { print } } }' jwt_2025.txt >> special.txt
  fi

# Add EUCLID TLE from 2025 list (see https://github.com/Bill-Gray/tles/ => 23092a25.tle)
  if [ ! -f $OUTDIR/euclid_2025.txt ]; then
	echo To have JWT TLEs move euclid_2025.txt from the scripts to the TLE dir.
  else
# Use approx int MJD. See also the C code in src directory.
	MJD=`expr $(date +%s) / 86400 + 40588`
	awk -v mjd="${MJD}" '$0 ~ mjd { print ; for(n=0; n<3; n++) { getline ; if ( index($1, "Euclid") > 0 ) { print "EUCLID"} else { print } } }' euclid_2025.txt >> special.txt
  fi

# Add JUICE TLE from 2024 list (see https://github.com/Bill-Gray/tles/ => 23053a.tle - 2025 not present)
  if [ ! -f $OUTDIR/juice_2024.txt ]; then
	echo To have JUICE TLEs move juice_2024.txt from the scripts to the TLE dir.
  else
# Use approx int MJD. See also the C code in src directory.
	MJD=`expr $(date +%s) / 86400 + 40588`
	awk -v mjd="${MJD}" '$0 ~ mjd { print ; for(n=0; n<3; n++) { getline ; if ( index($1, "JUICE") > 0 ) { print "JUICE"} else { print } } }' juice_2024.txt >> special.txt
  fi


#
# Produce the updated list of number of sats in the TLE file. Comment out if not needed.
#
$BINDIR/tle_satcount.sh

date | tee -a $LOGFILE
