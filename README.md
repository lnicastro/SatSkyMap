# SatSkyMap
A tool to find artificial satellites within a given sky region (RA/Dec - J2000),
as seen from a given point on Earth at a given epoch.
Single satellite info, regardless of its sky position, can be requested, too.
The executable tool name is `sat_skymap`.

[![DOI](https://zenodo.org/badge/256285121.svg)](https://zenodo.org/badge/latestdoi/256285121)

An interactive [web tool](https://sats.oas.inaf.it/) is available at INAF-OAS Bologna. This site also privides direct HTTPS queries via PHP wrapper.

An additional tool to check for a given satellite visibility is provided too: `sat_visibility`. It can be invoked via HTTPS too.
See below.

## Requirements
Two Lines Element (TLE) files, possibly with *up-to-date parameters*. See [here](https://celestrak.com/columns/v04n03/) for a description.

Two files, `default.tle` and `stations.txt`, are provided for testing purpose only.

The shell script `tle_retrieve.sh` (in `scripts`) is provided for automatic retrieval of TLE files from [CeleStrak](https://celestrak.com/).
Edit it to match your needs.
The two files `ALL_merged.lis` and `ALL_merged_nodeb.lis` list the retrieved files, including/excluding the debris objects, respectively.

## Compile and install
Code meant for Unix/Linux OS (including Mac OS). Can be easily adapted for Windows.
Just run `make` and move `sat_skymap` to your preferred folder. 

## Procedure
Use the Near-Earth type (SGP4) and Deep-Space type Ephemeris (SDP4) satellite motion model
to find those present in a given sky region at a given time observing from a given site on Earth.

`sat_skymap` requires as input:
1. a TLE file (*first or last command-line argument*),
2. the observer position (geodetic Lon, Lat, Alt),
3. the (circular) region center (RA, Dec - J2000 - can be omitted if geodetic position used as reference),
4. the region search radius,
5. reference date-time.

Invoke the tool without parameters to see predefined parameters and example. Use the `-h` switch for help:

```shell
% ./sat_skymap -h
Usage:
  sat_skymap tle_file [OPTIONS]

OPTIONS are:
  -a Alt_min,Alt_max  Geodetic altitude range filter (km; -G assumed by def.)
  -d CalendarDate     Calendar date (UTC) of interest (in the form yyyy-mm-ddThh:mm:ss[.sss])
  -D DirTLEs          Directory with the repository of TLE files (def. ./; ignore -T option)
  -i SatIntnlName     Single satellite selection via its international designator (region ignored)
  -j MJD              Modified Julian Date of interest (ignored if Calendar Date given)
  -l Lat,Lon,Alt      Geodetic observing site (comma separated data with no spaces)
  -n MaxSats          Maximum number of satellites to return (def. 1000)
  -p RA,Dec           J2000 sky coordinates of region to check
  -r radius           Region radius centered at the given coords
  -s SatNorad_n       Single satellite selection via its NORAD number (region ignored)
  -S SatName          Satellites selection via string name (substring matching applies; region ignored)
  -t deltaT           Second epoch delta time (seconds; def. 1)

Switches:
  -h            print this help
  -E            Use input geodetic location as reference location for region (-G by def.; -p ignored)
  -G            Compute (and return) geodetic location for each satellite
  -H            Compute sky separation via Haversine formula rather than cartesian triangle (suggested!)
  -I            Only information about the returned data and number of satellites found
                ('satellites' object not returned)
  -T            Use default repository directory for TLE files (def. ./; see sat_skymap_def.h) 
  -V            Return data only for sunlit satellites (potentially visible)

Example usage:
  sat_skymap default.tle -l-29.25627,-70.73805,2400 -p90.5,-30.3 -j58861.5 -r20
  sat_skymap stations.txt -H -l44.52804,11.33715,23.5 -j58941.71931 -s 25544
```

The space after the option character is optional.

Can give a file with a "list" of TLE files adding a `@` before the file name, e.g.:
```
  sat_skymap @ALLsats.lis -H -D/usr/local/TLErepo -l44.52804,11.33715,23.5 -j58980.5 -s 25989
```

## Output
A JSON string with information about the satellites in the FoV (see below).

## Examples
**ESO La Silla Observatory**

Scan the TLE element file 'default.tle' for satellites visible from
-  (latitude, longitude, altitude) = (-29.25627, -70.7380, 2400),
-  on MJD 58861.5 (UTC: 2020-01-13 12h = 2020-01-13T12:00:00),
-  at (RA, Dec) = 90.5, +30.3 (deg), within a 20-deg search radius.

Positions at given JD + deltaT (def. 1 s) is also computed and returned.
Additional computed info:
-  Local Mean Sidereal Time and Greenwich Mean Sidereal Time (hours)
-  Region Az, Alt, Parang
-  Sun RA/Dec, Alt/Az, geodetic longitude coordinates, parallactic angle and distance from region center (or zenith).
-  Geodetic Lon, Lat, Alt and theta (*for single satellite enquire* or if requested)

```
./sat_skymap default.tle -l-29.25627,-70.73805,2400 -j58861.5 -p90.5,-30.3 -r20
./sat_skymap default.tle -l-29.25627,-70.73805,2400 -d2020-01-13T12:00:00 -p90.5,-30.3 -r20

  {
  "swinfo": {"name": "sat_skymap", "author": "L. Nicastro @ INAF-OAS", "date": "2024-04-22", "version": "0.3f"},
  "input_params": {"tle_file": "default.tle", "location": ["lat":-29.2563, "lon": -70.7381, "alt":  2400.0],
    "region": {"ra":  90.5000, "dec":-30.3000, "radius": 20.0000, "lmst": 14.7803, "az": 222.1310, "alt":-14.4561, "parang": 137.324},
    "mjd": 58861.50000, "epoch_UTC": "2020-01-13T12:00:00", "gmst": 19.4962, "delta_time_s": 1, "max_sats": 1000,
    "notes": "All coordinates and radius in degrees. GMST, LMST in hrs."},
  "sun": {"ra":294.566, "dec":-21.517, "az": 101.818, "alt": 24.735, "lon":   2.123, "parang":-113.376, "separation_deg":123.255},
  "data_fields": {"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "in_sunlit", "HPXID_8"],
    "desc": ["RA Tinit", "Dec Tinit", "RA Tend", "Dec Tend", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "sunlit sat. flag", "HEALPix order 8 nested schema ID"],
    "type": ["double", "double", "double", "double", "double", "float", "float", "float", "int", "int"],
    "unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", ""]},
  "satellites": [{"name": "STARLINK-23", "intl_desig": "1919-029C ", "norad_n": 44237,
    "data": [ 77.5770, -17.4894,  77.6262, -17.4900,   7660.9,18.1955,  90,  2.818,1,339078]},

  { ... }

  }],
  "run_command": "sat_skymap default.tle -l-29.25627,-70.73805,2400 -p90.5,-30.3 -j58861.5 -r20", "status": 0, "errmsg": "", "n_sats_found": 6, "n_sats": 6, "n_sunlit_sats": 6
  }
```
It is:
-  Distance: distance to satellite in *km*
-  Separation: angular separation from the search point in *degrees*
-  PA: position angle of motion in *degrees*
-  Speed: apparent angular rate of motion in *arcminutes/second* (== *degrees/minute*)
-  sunlit: *0* => in Earth shade, *1* => sunlit (preliminary; work in progress)
-  HPXID_8: HEALPix order 8 nested schema ID


**ISS position (NORAD ID = 25544)**

JSON shown in expanded format.
```
./sat_skymap stations.txt -H -l44.52804,11.33715,23.5 -j58959.53 -s 25544

{
	"swinfo": {
		"name": "sat_skymap",
		"author": "L. Nicastro @ INAF-OAS",
		"date": "2024-04-22",
		"version": "0.3f"
	},
	"geoloc_fields":{
		"lat":{
			"desc":"Geodetic Latitude",
			"unit":"deg"
		},
		"lon":{
			"desc":"Geodetic Longitude",
			"unit":"deg"
		},
		"alt":{
			"desc":"Geodetic Altitude",
			"unit":"km"
		},
		"theta":{
			"desc":"Equatorial angle (Lon + GMST = RA)",
			"unit":"deg"
		}
	},
	"input_params": {
		"tle_file": "stations.txt",
		"location": {
			"lat": 44.5280,
			"lon": 11.3371,
			"alt": 23.5
		},
		"region": {
			"ra": 51.2026,
			"dec": 44.5280,
			"radius": 20.0000,
			"lmst": 3.4135,
			"az": 0.0000,
			"alt": 90.0000,
			"parang": 180.000
		},
		"mjd": 58959.53000,
		"epoch_UTC": "2020-04-20T12:43:12",
		"gmst": 2.6577,
		"delta_time_s": 1,
		"max_sats": 1000,
		"notes": "All coordinates and radius in degrees. GMST, LMST in hrs."
	},
	"sun": {
		"ra": 28.769,
		"dec": 11.785,
		"az": 217.381,
		"alt": 52.026,
		"lon": -11.096,
		"parang": 26.240,
		"separation_deg": 37.974
	},
	"data_fields": {
		"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "in_sunlit", "HPXID_8"],
		"desc": ["RA T_ini", "Dec T_ini", "RA T_end", "Dec T_end", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "sunlit sat. flag", "HEALPix order 8 nested schema ID"],
		"type": ["double", "double", "double", "double", "double", "float", "float", "float", "int", "int"],
		"unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", "", ""]
	},
	"satellites": [{
		"name": "ISS (ZARYA)",
		"intl_desig": "1998-067A ",
		"norad_n": 25544,
		"geoloc": {
			"lat": -44.174,
			"lon": 103.841,
			"alt": 433.24,
			"theta": 143.706
		},
		"data": [185.2172, -53.2256, 185.2753, -53.2252, 11436.45, 149.1219, 89, 2.087, 0, 690893]
	}],
	"run_command": "sat_skymap stations.txt -H -l44.52804,11.33715,23.5 -j58959.53 -s 25544",
	"status": 0,
	"errmsg": "",
	"n_sats_found": 1,
	"n_sats": 1
	"n_sunlit_sats": 0
}
```
**Search satellite(s) by (initial) name**

JSON shown in expanded format.

```
./sat_skymap stations.txt -S LEMUR

{
   "swinfo": {
	"name": "sat_skymap",
	"author": "L. Nicastro @ INAF-OAS",
	"date": "2024-04-22",
	"version": "0.3f"
   },
   "input_params": {
	"tle_file": "stations.txt",
	"location":{
         "lat": -29.2563,
         "lon": -70.7380,
         "alt": 2400.0
	},
	"region": {
         "ra": 90.5000,
         "dec": -30.3000,
         "radius": 20.0000,
         "lmst": 14.7803,
         "az": 222.1309,
         "alt": -14.4561,
         "parang": 137.324
	},
	"mjd": 58861.50000,
	"epoch_UTC": "2020-01-13T12:00:00",
	"gmst": 19.4962,
	"delta_time_s": 1,
	"max_sats": 1000,
	"notes": "All coordinates and radius in degrees. GMST, LMST in hrs."
   },
   "sun": {
	"ra": 294.566,
	"dec": -21.517,
	"az": 101.818,
	"alt": 24.735,
	"lon": 2.123,
	"parang": -113.376,
	"separation_deg": 123.255
   },
   "data_fields": {
	"name": [
		"RA_start",
		"Dec_start",
		"RA_end",
		"Dec_end",
		"Distance",
		"Separation",
		"PA",
		"Speed",
		"in_sunlit",
		"HPXID_8"
	],
	"desc": [
		"RA T_ini",
		"Dec T_ini",
		"RA T_end",
		"Dec T_end",
		"distance to sat.",
		"angular separation",
		"position angle",
		"apparent angular rate of motion",
		"sunlit sat. flag",
		"HEALPix order 8 nested schema ID"
	],
	"type": [
		"double",
		"double",
		"double",
		"double",
		"double",
		"float",
		"float",
		"float",
		"int"
		"int"
	],
	"unit": [
		"deg",
		"deg",
		"deg",
		"deg",
		"km",
		"deg",
		"deg",
		"arcmin/s",
		""
		""
	]
   },
   "satellites": [
 	{
		"name": "LEMUR-2-VU",
		"intl_desig": "2018-046E",
		"norad_n": 43558,
		"data": [ 90.0544, 75.2547, 89.9533, 75.2957, 8573.77, 105.5557, 327, 2.904, 1, 126719 ]
	},
	{
		"name": "LEMUR-2-ALEXANDER",
		"intl_desig": "2018-046F",
		"norad_n": 43559,
		"data": [ 351.8948, 17.6476, 351.9209, 17.6085, 9443.62, 265.7559, 147, 2.781, 1, 321031 ]
	},
	{
		"name": "LEMUR-2-YUASA",
		"intl_desig": "2018-046G",
		"norad_n": 43560,
		"data": [ 85.0062, 8.0375, 85.0374, 8.0723, 9532.90, 38.7291, 41, 2.792, 0, 379090 ]
	},
	{
		"name": "LEMUR-2-TOMHENDERSON",
		"intl_desig": "2018-046H",
		"norad_n": 43561,
		"data": [ 77.0633, 0.0105, 77.0993, 0.0394, 9723.66, 33.1553, 51, 2.772, 0, 366949 ]
	}
   ],
   "run_command": "sat_skymap stations.txt -S LEMUR",
   "status": 0,
   "errmsg": "",
   "n_sats_found": 4,
   "n_sats": 4
   "n_sunlit_sats": 2
}
```

An example of direct HTTPS query:
```
https://sats.oas.inaf.it/get_sats.php?tlename=active&location=44.523,11.339,0&epoch=2024-08-09T12:47:33&deltat=10&coords=161.733,44.523&size_x=19.840&maxsats=2000&geoc=1
```

## Satellite visibility
The tool `sat_visibility` can be used to check for a satellite visibility (sunlit) at a given site given its NORAD number or international designator.
```shell
% ./sat_visisbility -h
Usage:
  sat_visibility tle_file [OPTIONS]

OPTIONS are:
  -a Alt_min,Alt_max Geodetic altitude range filter (km)
  -d CalendarDates   Calendar date (UTC) range of interest (in the form yyyy-mm-ddThh:mm:ss[.sss],yyyy-...)
  -D DirTLEs         Directory with the repository of TLE files (def. ./; ignore -T option)
  -i SatIntnlName    Satellite selection via its international designator
  -j MJDs            Modified Julian Date range of interest (MJDstart,MJDend - ignored if Calendar Date given)
  -l Lat,Lon,Alt     Geodetic observing site (comma separated data with no spaces)
  -r radius          Region search radius - centred at the site zenith (degrees; def. 90)
  -s SatNorad_n      Satellite selection via its NORAD number
  -t deltaT          Visibility check delta time (seconds; def. 3)

Switches:
  -h                 Print this help
  -N                 Skip visibility check until Sun is set - if it is up at start date
  -T                 Use default repository directory for TLE files (def. ./; see sat_visibility_def.h)

Example usage:
  sat_visibility stations.txt -T -l44.52324,11.33914,23 -d2024-05-11T20:00:00,2024-05-11T22:00:00 -s 25544
  sat_visibility stations.txt -T -t 10 -l-29.25627,-70.73805,2400 -j60436.54,60436.55 -s 25544
```

A JSON output string gives all the relevant information about the satellite and the sun position, in particucular:
```
    Times of Start/End/Mid of the visivility range.
    The highest satellite altitude over the horizon.

    Distance to satellite in km
    Angular separation from the zenith
    Position angle of motion
    Apparent angular rate of motion in arcminutes/second
```

An example of direct HTTPS query:
```
https://sats.oas.inaf.it/get_visibility.php?tlename=active&location=44.523,11.339,0&epoch=2024-08-09T12:47:33&norad_n=43571&view
```

See the [web tool](https://sats.oas.inaf.it/) for automatically created query links.


## Acknowledgements
Parts of the code from [Project Pluto](https://www.projectpluto.com/), Bill Gray's `sat_code` on [GitHub](https://github.com/Bill-Gray/sat_code).

TLE files from [CeleStrak](https://celestrak.com/) and [Space-Track](https://www.space-track.org/).

HEALPix C library used. See [HEALPix](https://healpix.sourceforge.io/).
