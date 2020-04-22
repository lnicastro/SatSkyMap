# SatSkyMap
A tool to find artificial satellites within a given sky region (RA/Dec - J2000),
as seen from a given point on Earth at a given epoch.
Single satellite info, regardless of its position can be requested, too.
The executable tool name is `sat_skymap`.

An interactive [web tool](https://cats.oas.inaf.it/SatSkyweb/) is available at INAF OAS Bologna.

## Requirements
Two Lines Element (TLE) files, possibly with *up-to-date parameters*. See [here](https://celestrak.com/columns/v04n03/) for a description.

Two files, `default.tle` and `stations.txt`, are provided for testing purpose only.

The shell script `tle_retrieve.sh` (in `scripts`) is provided for automatic retrieval of TLE files from [CeleStrak](https://celestrak.com/).
The two files `ALL_merged.lis` and `ALL_merged_nodeb.lis` list the retrieved files, including/excluding the debris objects, respectively.

## Compile and install
Code meant for Unix/Linux OS (including Mac OS). Can be easily adapted for Windows.
Just run `make` and move `sat_skymap` to your preferred folder. 

## Procedure
Use the Near-Earth type (SGP4) and Deep-Space type Ephemeris (SDP4) satellite motion model
to find those present in a given sky region at a given time observing from a given site on Earth.

`sat_skymap` requires as input:
1. a TLE file (*first command-line argument*),
2. the observer position (geodetic Lon, Lat, Alt),
3. the (circular) region RA/Dec (J2000) and search radius,
4. date/time.

Invoke the tool without parameters to see predefined parameters and example. Use the `-h` switch for help:

```shell
% ./sat_skymap -h
Usage:
  sat_skymap tle_file [OPTIONS]

OPTIONS are:
  -d CalendarDate	Calendar date (UTC) of interest (in the form yyyy-mm-ddThh:mm:ss[.sss])
  -i sat_intnlname	(TODO) Single satellite selection via its international designator (region ignored)
  -j MJD		Modified Julian Date of interest (ignored if Calendar Date given)
  -l Lat,Lon,Height	Observation site (comma separated data with no spaces)
  -n MaxSats		Maximum number of satellites to return (def. 1000)
  -p RA,Dec		J2000 sky coordinates of region to check
  -r radius		Region radius centered at the given coords
  -s sat_norad_n	Single satellite selection via its NORAD number (region ignored)
  -t deltat		Second epoch delta time (s, def. 1)
  -D DirTLEs		Directory with the repository of TLE files (def ./; ignore -T option)

Switches:
  -h			print this help
  -H			Compute sky separation via Haversine formula rather than cartesian triangle (suggested! Default?)
  -I			Information about the returned data and number of satellites found
 			('satellites' object not returned)
  -T			Use default repository directory for TLE files (see sat_skymap_def.h; def. ./) 

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
-  (lat, lon, height) = (-29.25627, -70.7380, 2400),
-  on MJD 58861.5 (UTC: 2020-01-13 12h = 2020-01-13T12:00:00),
-  at (RA, Dec) = 90.5, +30.3 (deg), within a 20-deg search radius.

Positions at given JD + deltat (def. 1) s is also computed and returned.
Additional computed info:
-  Local Mean Sidereal Time and Greenwich Mean Sidereal Time (hours)
-  Region Az, Alt, Parang
-  Sun RA/Dec, Alt/Az coordinates
-  Geodetic Lon, Lat, Alt and theta (*for single satellite enquire*)

```
./sat_skymap default.tle -l-29.25627,-70.73805,2400 -j58861.5 -p90.5,-30.3 -r20
./sat_skymap default.tle -l-29.25627,-70.73805,2400 -d2020-01-13T12:00:00 -p90.5,-30.3 -r20

  {
  "swinfo": {"name": "sat_skymap", "author": "L. Nicastro @ INAF-OAS", "date": "2020-04-20", "version": "0.2b"},
  "input_params": {"tle_file": "default.tle", "location": ["lat":-29.2563, "long": -70.7381, "height":  2400.0],
    "region": {"ra":  90.5000, "dec":-30.3000, "radius": 20.0000, "lmst": 14.7803, "az": 222.1310, "alt":-14.4561, "parang": 137.324},
    "mjd": 58861.50000, "epoch_UTC": "2020-01-13T12:00:00", "gmst": 19.4962, "delta_time_s": 1, "max_sats": 1000,
    "notes": "All coordinates and radius in degrees. GMST, LMST in hrs."},
  "sun": {"ra":294.566, "dec":-21.517, "az": 101.818, "alt": 24.735, "parang":-113.376, "separation_deg":123.255},
  "data_fields": {"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "HPXID_8"],
    "desc": ["RA Tinit", "Dec Tinit", "RA Tend", "Dec Tend", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "HEALPix order 8 nested schema ID"],
    "type": ["double", "double", "double", "double", "double", "float", "float", "float", "int"],
    "unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", ""]},
  "satellites": [{"name": "STARLINK-23", "intl_desig": "1919-029C ", "norad_n": 44237,
    "data": [ 77.5770, -17.4894,  77.6262, -17.4900,   7660.9,18.1955,  90,  2.818,339078]},

  { ... }

  }],
  "status": 0, "errmsg": "", "n_sats_found": 6, "n_sats": 6
  }
```
It is:
-  Distance: distance to satellite in *km*,
-  Separation: angular separation from the search point in *degrees*,
-  PA: position angle of motion in *degrees*,
-  Speed: apparent angular rate of motion in *arcminutes/second* (== *degrees/minute*),
-  HPXID_8: HEALPix order 8 nested schema ID


**ISS position (NORAD ID = 25544)**

JSON shown in expanded format.
```
./sat_skymap stations.txt -H -l44.52804,11.33715,23.5 -j58959.53 -s 25544

{   "swinfo":{
      "name":"sat_skymap",
      "author":"L. Nicastro @ INAF-OAS",
      "date":"2020-04-20",
      "version":"0.2b"
   
},
   "input_params":{  "tle_file":"stations.txt",
      "location":{
         "lat":44.5280,
         "long":11.3371,
         "height":23.5
      
},
      "region":{
         "ra":51.2026,
         "dec":44.5280,
         "radius":20.0000,
         "lmst":3.4135,
         "az":0.0000,
         "alt":90.0000,
         "parang":180.000
      
},
      "mjd":58959.53000,
      "epoch_UTC":"2020-04-20T12:43:12",
      "gmst":2.6577,
      "delta_time_s":1,
      "max_sats":1000,
      "notes":"All coordinates and radius in degrees. GMST, LMST in hrs."
   
},
   "sun":{
      "ra":28.769,
      "dec":11.785,
      "az":217.381,
      "alt":52.026,
      "parang":26.240,
      "separation_deg":37.974
   
},
   "geoloc":{
      "lat":-44.174,
      "lon":103.841,
      "alt":433.24,
      "theta":143.706
   
},
   "data_fields":{  "name":[
         "RA_start",
         "Dec_start",
         "RA_end",
         "Dec_end",
         "Distance",
         "Separation",
         "PA",
         "Speed",
         "HPXID_8"
      
	],
      "desc":[
         "RA T_ini",
         "Dec T_ini",
         "RA T_end",
         "Dec T_end",
         "distance to sat.",
         "angular separation",
         "position angle",
         "apparent angular rate of motion",
         "HEALPix order 8 nested schema ID"
      
	],
      "type":[
         "double",
         "double",
         "double",
         "double",
         "double",
         "float",
         "float",
         "float",
         "int"
      
	],
      "unit":[
         "deg",
         "deg",
         "deg",
         "deg",
         "km",
         "deg",
         "deg",
         "arcmin/s",
         ""
	]
   
},
   "satellites":[  {  "name":"ISS (ZARYA)",
         "intl_desig":"1998-067A ",
         "norad_n":25544,
         "data":[
            185.2172,
            -53.2256,
            185.2753,
            -53.2252,
            11436.45,
            149.1219,
            89,
            2.087,
            690893
	]
      
}
   
],
   "status":0,
   "errmsg":"",
   "n_sats_found":1,
   "n_sats":1
}
```

## Acknowledgements
Parts of the code from [Project Pluto](https://www.projectpluto.com/), Bill Gray's `sat_code` on [GitHub](https://github.com/Bill-Gray/sat_code).

TLE files from [CeleStrak](https://celestrak.com/).

HEALPix C library used. See [HEALPix](https://healpix.sourceforge.io/).
