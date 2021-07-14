/*
   sat_skymap.c

  Purpose:
    Find which satellites are within a given angular separation (approx)
    of a given RA/Dec (J2000), as seen from a given point.
    Single satelletite info, regardless of its position can be requested, too.

  Procedure:
    The code reads in a TLE file (name provided as the first or last command-line
    argument) with details of the observer position, search radius, date/time,
    and RA/Dec provided on the command line.

  Input Parameters:
    -a Alt_min,Alt_max	Geodetic altitude range filter (km; -G assumed by def.)
    -d CalendarDate	Calendar date (UTC) of interest (in the form yyyy-mm-ddThh:mm:ss[.sss])
    -D DirTLEs		Directory with the repository of TLE files (def. ./; ignore -T option)
    -i SatIntnlName	Single satellite selection via its international designator (region ignored)
    -j MJD		Modified Julian Date of interest (ignored if Calendar Date given)
    -l Lat,Lon,Alt	Geodetic observing site (comma separated data with no spaces)
    -n MaxSats		Maximum number of satellites to return (def. 1000)
    -p RA,Dec		J2000 sky coordinates of region to check
    -r radius		Region radius centered at the given coords
    -s SatNorad_n	Single satellite selection via its NORAD number (region ignored)
    -S SatName		Satellites selection via string name (substring matching applies; region ignored)
    -t deltaT		Second epoch delta time (seconds; def. 1)

  Switches:
    -h			print this help
    -E			Use input geodetic location as reference location for region (-G by def.; -p ignored)
    -G			Compute (and return) geodetic location for each satellite
    -H			Compute sky separation via Haversine formula rather than cartesian triangle (suggested! Default?)
    -I			Only information about the returned data and number of satellites found
  			('satellites' object not returned)
    -T			Use default repository directory for TLE files (def. ./; see sat_skymap_def.h) 
    -V			Return data only for satellites in sunlight (potentially visible)

  Output:
    A json string with information about the satellites in the FoV (see below).

  Example (La Silla coords):

    ./sat_skymap default.tle -l-29.25627,-70.73805,2400 -j58861.5 -p90.5,-30.3 -r20
    ./sat_skymap default.tle -l-29.25627,-70.73805,2400 -d2020-01-13T12:00:00 -p90.5,-30.3 -r20

  Scan the TLE element file 'default.tle' for satellites visible from
  (latitude, longitude, altitude) = (-29.25627, -70.7380, 2400),
  on MJD 58861.5 (UTC: 2020-01-13 12h = 2020-01-13T12:00:00),
  at (RA, Dec) = 90.5, +30.3 (deg), within a 20-deg search radius.
  Positions at given JD + deltaT (def. 1 s) is also computed and returned.
  Additional computed info:
    Local Mean Sidereal Time
    region Az, Alt, Parang
    Sun RA/Dec, Alt/Az, geodetic longitude coordinates, parallactic angle
    and distance from region center (or zenith).

  The output looks like this:

  {
  "swinfo": {"name": "sat_skymap", "author": "L. Nicastro @ INAF-OAS", "date": "2021-07-06", "version": "0.3b"},
  "input_params": {"tle_file": "default.tle", "location": ["lat":-29.2563, "lon": -70.7381, "alt":  2400.0],
    "region": {"ra":  90.5000, "dec":-30.3000, "radius": 20.0000, "lmst": 14.7803, "az": 222.1310, "alt":-14.4561, "parang": 137.324},
    "mjd": 58861.50000, "epoch_UTC": "2020-01-13T12:00:00", "gmst": 19.4962, "delta_time_s": 1, "max_sats": 1000,
    "notes": "All coordinates and radius in degrees. GMST, LMST in hrs."},
  "sun": {"ra":294.566, "dec":-21.517, "az": 101.818, "alt": 24.735, "lon":   2.123, "parang":-113.376, "separation_deg":123.255},
  "data_fields": {"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "in_sunlight", "HPXID_8"],
    "desc": ["RA Tinit", "Dec Tinit", "RA Tend", "Dec Tend", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "sat. in sunlight flag", "HEALPix order 8 nested schema ID"],
    "type": ["double", "double", "double", "double", "double", "float", "float", "float", "int", "int"],
    "unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", "", ""]},
  "satellites": [{"name": "STARLINK-23", "intl_desig": "1919-029C ", "norad_n": 44237,
    "data": [ 77.5770, -17.4894,  77.6262, -17.4900,   7660.9,18.1955,  90,  2.818,1,339078]},

  { ... }

  }],
  "status": 0, "errmsg": "", "n_sats_found": 6, "n_sats": 6, "n_sats_in_sunlight": 6
  }

  Where:
    Distance: distance to satellite in km
    Separation: angular separation in degrees from the search point
    PA: position angle of motion
    Speed: apparent angular rate of motion in arcminutes/second (or degrees/minute)
    in_sunlight: 0 => in Earth shade, 1 => in sunlight (preliminary; work in progress)
    HPXID_8: HEALPix order 8 nested schema ID

  Example (ISS position):
    ./sat_skymap stations.txt -H -l44.52804,11.33715,23.5 -j58959.53 -s 25544

  Output:

  {
  "swinfo": {"name": "sat_skymap", "author": "L. Nicastro @ INAF-OAS", "date": "2021-07-06", "version": "0.3b"},
  "geoloc_fields": {"lat": {"desc": "Geodetic Latitude", "unit": "deg"}, "lon": {"desc": "Geodetic Longitude", "unit": "deg"}, "alt": {"desc": "Geodetic Altitude", "unit": "km"}, "theta": {"desc": "Equatorial angle (Lon + GMST = RA)", "unit": "deg"}},
  "input_params": {"tle_file": "stations.txt", "location": {"lat": 44.5280, "lon":  11.3371, "alt":    23.5},
    "region": {"ra":  51.2026, "dec": 44.5280, "radius": 20.0000, "lmst":  3.4135, "az":   0.0000, "alt": 90.0000, "parang": 180.000},
    "mjd": 58959.53000, "epoch_UTC": "2020-04-20T12:43:12", "gmst":  2.6577, "delta_time_s": 1, "max_sats": 1000,
    "notes": "All coordinates and radius in degrees. GMST, LMST in hrs."},
  "sun": {"ra": 28.769, "dec": 11.785, "az": 217.381, "alt": 52.026, "lon": -11.096, "parang":  26.240, "separation_deg": 37.974},
  "data_fields": {"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "HPXID_8"],
  "desc": ["RA T_ini", "Dec T_ini", "RA T_end", "Dec T_end", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "HEALPix order 8 nested schema ID"],
  "type": ["double", "double", "double", "double", "double", "float", "float", "float", "int"],
  "unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", ""]},
  "satellites": [{"name": "ISS (ZARYA)", "intl_desig": "1998-067A ", "norad_n": 25544,
    "geoloc": {"lat":-44.174, "lon":103.841, "alt":   433.24, "theta": 143.706},
    "data": [185.2172,-53.2256,185.2753,-53.2252, 11436.45,149.1219, 89, 2.087,690893]}],
  "status": 0, "errmsg": "", "n_sats_found": 1, "n_sats": 1, "n_sats_in_sunlight": 0
  }


  LN @ INAF-OAS, Jan 2020.  Last change: 14/07/2021
*/

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "norad.h"
#include "observe.h"

#include "const_def.h"
#include "chealpix.h"
#include "sat_skymap_def.h"


/*
  From Spherical coords in radians (e.g. RA, Dec) to unit vector
*/

void sphrad2v(const double rarad, const double decrad, double v[3])
{
  double cosd = cos(decrad);

  v[0] = cos(rarad) * cosd;
  v[1] = sin(rarad) * cosd;
  v[2] = sin(decrad);
}


/*
  Trim trailing space
*/
char *trimend(char *str)
{
  char *end;

  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end)) end--;
  end[1] = '\0';

  return str;
}


void Usage() {
  printf("Usage:\n  %s tle_file [OPTIONS]\n\n"
  "OPTIONS are:\n"
  "  -a Alt_min,Alt_max	Geodetic altitude range filter (km; -G assumed by def.)\n"
  "  -d CalendarDate	Calendar date (UTC) of interest (in the form yyyy-mm-ddThh:mm:ss[.sss])\n"
  "  -D DirTLEs		Directory with the repository of TLE files (def. ./; ignore -T option)\n\n"
  "  -i SatIntnlName	Single satellite selection via its international designator (region ignored)\n"
  "  -j MJD		Modified Julian Date of interest (ignored if Calendar Date given)\n"
  "  -l Lat,Lon,Alt	Geodetic observing site (comma separated data with no spaces)\n"
  "  -n MaxSats		Maximum number of satellites to return (def. 1000)\n"
  "  -p RA,Dec		J2000 sky coordinates of region to check\n"
  "  -r radius		Region radius centered at the given coords\n"
  "  -s SatNorad_n	Single satellite selection via its NORAD number (region ignored)\n"
  "  -S SatName		Satellites selection via string name (substring matching applies; region ignored)\n"
  "  -t deltaT		Second epoch delta time (seconds; def. 1)\n"
  "\nSwitches:\n"
  "  -h			print this help\n"
  "  -E			Use input geodetic location as reference location for region (-G by def.; -p ignored) \n"
  "  -G			Compute (and return) geodetic location for each satellite\n"
  "  -H			Compute sky separation via Haversine formula rather than cartesian triangle (suggested!)\n"
  "  -I			Only information about the returned data and number of satellites found\n"
  " 			('satellites' object not returned)\n"
  "  -T			Use default repository directory for TLE files (def. ./; see sat_skymap_def.h)\n"
  "  -V			Return data only for satellites in sunlight (visible)\n\n"

  "Example usage:\n"
  "  %s default.tle -l-29.25627,-70.73805,2400 -p90.5,-30.3 -j58861.5 -r20\n"
  "  %s stations.txt -H -l44.52804,11.33715,23.5 -j58941.71931 -s 25544\n", progname, progname, progname);
}


void open_json() {
  printf("{"
    "\"swinfo\": {"
    "\"name\": \"%s\", \"author\": \"%s\", \"date\": \"%s\", \"version\": \"%s\"}, ",
    progname, progauthor, progdate, progversion);
}

void add_params_json(Params p) {
  printf("\"input_params\": {"
    "\"tle_file\": \"%s\", \"location\": {\"lat\":%8.4lf, \"lon\":%9.4lf, \"alt\":%8.1lf}, \"region\": {\"ra\":%9.4lf, \"dec\":%8.4lf, \"radius\":%8.4lf, \"lmst\":%8.4lf, \"az\":%9.4lf, \"alt\":%8.4lf, \"parang\":%8.3lf}, "
    "\"mjd\": %11.5lf, \"epoch_UTC\": \"%s\", \"gmst\":%8.4lf, \"delta_time_s\": %d, \"max_sats\": %d, "
    "\"notes\": \"All coordinates and radius in degrees. GMST, LMST in hrs.\"}, ",
    p.tle_file_name, p.lat, p.lon, p.ht_in_meters, p.ra_deg, p.de_deg, p.search_radius, p.lmst, p.az, p.alt, p.parang, p.mjd, p.date, p.gmst, p.delta_time, p.max_sats);
}

void add_sundata_json(Sun sun) {
  printf("\"sun\": {"
    "\"ra\":%7.3lf, \"dec\":%7.3lf, \"az\":%8.3lf, \"alt\":%7.3lf, \"lon\":%8.3lf, \"parang\":%8.3lf, \"separation_deg\":%7.3lf}, ",
    sun.ra, sun.dec, sun.az, sun.alt, sun.lon, sun.parang, sun.sep);
}

void add_fieldsdesc_json() {
  printf(
    "\"data_fields\": {\"name\": [\"RA_start\", \"Dec_start\", \"RA_end\", \"Dec_end\", \"Distance\", \"Separation\", \"PA\", \"Speed\", \"in_sunlight\", \"HPXID_8\"], "
    "\"desc\": [\"RA T_ini\", \"Dec T_ini\", \"RA T_end\", \"Dec T_end\", \"distance to sat.\", \"angular separation\", \"position angle\", \"apparent angular rate of motion\", \"sat. in sunlight flag\", \"HEALPix order 8 nested schema ID\"], "
    "\"type\": [\"double\", \"double\", \"double\", \"double\", \"double\", \"float\", \"float\", \"float\", \"int\", \"int\"], "
    "\"unit\": [\"deg\", \"deg\", \"deg\", \"deg\", \"km\", \"deg\", \"deg\", \"arcmin/s\", \"\", \"\"]},");
}

void add_geofieldsdesc_json() {
  printf(
    "\"geoloc_fields\": {"
	"\"lat\": {\"desc\": \"Geodetic Latitude\", \"unit\": \"deg\"}, "
	"\"lon\": {\"desc\": \"Geodetic Longitude\", \"unit\": \"deg\"}, "
	"\"alt\": {\"desc\": \"Geodetic Altitude\", \"unit\": \"km\"}, "
	"\"theta\": {\"desc\": \"Equatorial angle (Lon + GMST = RA)\", \"unit\": \"deg\"}"
    "}, ");
}

void add_satlatlon_json(Geoloc geo) {
    printf( "\"geoloc\": {\"lat\":%7.3lf, \"lon\":%7.3lf, \"alt\":%9.2lf, \"theta\":%8.3lf}, ", geo.lat, geo.lon, geo.alt, geo.theta);
}

void close_stat_json(int status, char *errmsg, int n_sats_found, int n_sats, int n_sats_in_sunlight) {
  printf(" \"status\": %d, \"errmsg\": \"%s\", \"n_sats_found\": %d, \"n_sats\": %d,"
	 " \"n_sats_in_sunlight\": %d}\n", status, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
}


/* GMST in hours. Reference:  The 1992 Astronomical Almanac, page B6. */
/* Unused.
static inline double myThetaG(double jd)
{
  const double omega_E = 1.00273790934;  // Earth rotations per sidereal day (non-constant)
  const double UT = fmod(jd + .5, 1.);
  double t_cen, GMST;

  t_cen = (jd - UT - JD2000) / 36525.;
  GMST = 24110.54841 + t_cen * (8640184.812866 + t_cen * (0.093104 - t_cen * 6.2E-6));

  GMST = fmod( GMST + SEC_IN_DAY * omega_E * UT, SEC_IN_DAY);
  if( GMST < 0.)
     GMST += SEC_IN_DAY;

  return(24. * GMST / SEC_IN_DAY);
  }
*/


/* Compute intesection points distance on Earth of the Satellite --> Sun connecting line.
   Set to 0, 0 (and 0 returned) if Earth not intercepted (satellite is illuminated).

   See http://paulbourke.net/geometry/circlesphere/index.html#linesphere
   and intersect_line_and_sphere in SkyField python package. 
*/
int intersect_satsun_sphere(double *satpos, double *sunpos, double *eray)
{
  int i;
  double center[3], endpoint[3], l_endpoint = 0., minus_b = 0., c = 0., discriminant;

  for (i = 0; i < 3; i++) {
	center[i] = -satpos[i];
	endpoint[i] = sunpos[i] - satpos[i];

	c += center[i] * center[i];
	l_endpoint += endpoint[i] * endpoint[i];
  }

  c -= ERAD * ERAD;
  l_endpoint = sqrt(l_endpoint);

  for (i = 0; i < 3; i++)
	minus_b += endpoint[i]/l_endpoint * center[i];

  minus_b *= 2;

  discriminant = minus_b * minus_b - 4 * c;

  if ( discriminant < 0. ) {
	eray[0] = 0;
	eray[1] = 0;
	return(0);
  }

  double dsqrt = sqrt(discriminant);

  eray[0] = (minus_b - dsqrt) / 2.;
  eray[1] = (minus_b + dsqrt) / 2.;

  if ( eray[1] < 0 )
	return(0);

  return(1);
}



/* Compute satellite geodetic Lon, Lat, Alt and theta (RA) position given its ECI position and GMST.
   Returned values are deg, deg, km, deg.
*/
void sat_geoLocation(double gmst, double *pos,  Geoloc *geo) {
  double r, e2, phi, sinphi, c;

  geo->theta = atan2(pos[1], pos[0]) * RAD2DEG;  /* degrees */
  geo->lon = (geo->theta - gmst * 15.);
  if (geo->lon > 180)
	geo->lon -= 360.;
  else if (geo->lon < -180.)
	geo->lon += 360.;

  r = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  e2 = F * (2. - F);
  geo->lat = atan2(pos[2], r);

  do {
	phi = geo->lat;
	sinphi = sin(phi);
	c = 1. / sqrt(1. - e2 * sinphi * sinphi);
	geo->lat = atan2(pos[2] + ERAD * c * e2 * sinphi, r);
  } while (fabs(geo->lat - phi) >= 1e-10);

  geo->alt = r/cos(geo->lat) - ERAD * c;  /* kilometers */

  geo->lat *= RAD2DEG;  /* degrees */
   if (geo->lat > 90.)
	geo->lat -= 360.;  /* could keep in -90, 90 range, but then need to manage geo->lon accordingly */
}


int main(int argc, char **argv)
{
  char errmsg[128] = {""}, tle_path[128] = {"."}, line1[100], line2[100],
	sat_name[25], tle_path_file[200], opt, *endptr;
  double jd, rho_sin_phi, rho_cos_phi, observer_loc[3], observer_loc2[3],
	 target_ra, target_dec, target_lon = 0, target_lat = 0;
  int i, par_pos, status = 0, n_sats_found = 0, n_sats = 0, n_sats_in_sunlight = 0,
	 len_intl_desig = 0, len_satname = 0;
  bool in_region = false;  /* used for single satellite request */


  if ( argc < 2 ) {
	open_json();
	sprintf(errmsg, "Example usage:\n %s default.tle -l-29.25627,-70.73805,2400 -p90.5,-30.3 -j58861.5 -r20", progname);
 	close_stat_json(1, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
	exit(0);
  }


  Params p;

/* Set default values */
  p.lat = -29.25627;
  p.lon = -70.7380;
  p.ht_in_meters = 2400.;
  p.mjd = 58861.5;        /* 13 Jan 2020 12h UT */
  strcpy(p.date, "2020-01-13T12:00:00");
  p.ra_deg = 90.5;        /* Default search coords in degrees */
  p.de_deg = -30.3;
  p.search_radius = 20.;  /* Default search radius in degrees */
  p.delta_time = 1;       /* Default time step (s) for second epoch pos. computation */
  p.max_sats = 1000;      /* Maximum number of satellites entries to retrieve */
  p.haversine = false;    /* If to use Haversine formula for angular distance computation */
  p.info_only = false;    /* If just return info about json data and number of satellites in the requested region */
  p.single_sat_i = false; /* If enquire for just 1 satellite using its international designator */
  p.single_sat_n = false; /* If enquire for just 1 satellite using its NORAD number */
  p.satname_filter = false; /* If select via name matching string */
  p.geoloc_requested = false;  /* If geodetic location requested for each satellite */ 
  p.geoloc_reference = false;  /* If region selection must refer to geodetic location
				rather than RA, Dec. If true, p.geoloc_requested is also set to true. */ 
  p.use_deftledir = false;  /* If default directory with TLE files should be used (def. ./) */
  p.altrng_requested = false;  /* If altitude range filter requested */
  p.in_sunlight_only = false;  /* If only data for illuminated satellites requested */
  p.intl_desig[0] = '\0';
  p.satname[0] = '\0';
  p.norad_n = 0;
  p.alt_min = 0;
  p.alt_max = 1e32;


/* Note: no check on validity of input parameters! Could also use getarg. */

  for ( i = 1; i < argc; i++ )
    if ( argv[i][0] == '-' ) {
		opt = argv[i][1];
		par_pos = 2;
		if ( argv[i][2] == 0 ) {  /* blank after option */
//if ( argv[i+1][0] != '-' )  /* No, because of possible negative numbers after the param */
		  ++i;
		  par_pos = 0;
		}
	switch( opt )
	{
	  /* Switches */
 	  case 'h':
		Usage();
		exit(0);
 	  case 'E':
		p.geoloc_reference = true;
		p.geoloc_requested = true;
		i--;
		break;
 	  case 'G':
		p.geoloc_requested = true;
		i--;
		break;
 	  case 'H':
		p.haversine = true;
		i--;
		break;
 	  case 'I':
		p.info_only = true;
		i--;
		break;
 	  case 'T':
		p.use_deftledir = true;
		i--;
		break;
 	  case 'V':
		p.in_sunlight_only = true;
		i--;
		break;

	  /* Parameters */
          case 'a':
		sscanf(argv[i] + par_pos, "%lf,%lf", &p.alt_min, &p.alt_max);
		if ( p.alt_max < p.alt_min )  {
		  double dummy = p.alt_max;
		  p.alt_max = p.alt_min;
		  p.alt_min = dummy;
		}
		p.altrng_requested = true;
		p.geoloc_requested = true;
		break;
          case 'd':
		strcpy(p.date, argv[i] + par_pos);
		date2mjd(argv[i] + par_pos, &p.mjd);
		break;
 	  case 'D':
		strcpy(tle_path, argv[i] + par_pos);
		p.use_deftledir = false;
		break;
          case 'i':
		strcpy(p.intl_desig, argv[i] + par_pos);
		len_intl_desig = strlen(p.intl_desig);
		p.single_sat_i = true;
		p.altrng_requested = false;  // No altitude filter
		in_region = false;
		break;
          case 'j':
		p.mjd = atof(argv[i] + par_pos);
		strcpy(p.date, mjd2date(p.mjd));
		break;
          case 'l':
		sscanf(argv[i] + par_pos, "%lf,%lf,%lf", &p.lat, &p.lon, &p.ht_in_meters);
		/* Longitude locally in the range [0, 360[ (not [-180, +180]) degrees */
		if ( p.lon < 0. )
		  p.lon += 360;
		if ( p.lat < -90. || p.lat > 90. )  {
		  open_json();
		  sprintf(errmsg, "Latitude must be in the range [-90, +90] degrees. Read '%lf'", p.lat);
  		  close_stat_json(-1, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		  exit(-1);
		}
		break;
          case 'n':
		sscanf(argv[i] + par_pos, "%d", &p.max_sats);
		break;
          case 'p':
		sscanf(argv[i] + par_pos, "%lf,%lf", &p.ra_deg, &p.de_deg);
		if ( p.ra_deg < 0. || p.ra_deg > 360. )  {
		  open_json();
		  sprintf(errmsg, "RA must be in the range [0, 360[ degrees. Read '%lf'", p.ra_deg);
  		  close_stat_json(-1, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		  exit(-1);
		}
		if ( p.de_deg < -90. || p.de_deg > 90. )  {
		  open_json();
		  sprintf(errmsg, "Dec must be in the range [-90, +90] degrees. Read '%lf'", p.de_deg);
  		  close_stat_json(-1, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		  exit(-1);
		}
		in_region = true;
		break;
          case 'r':
		p.search_radius = atof(argv[i] + par_pos);
		break;
          case 's':
		p.norad_n = strtol(argv[i] + par_pos, &endptr, 10);
		p.single_sat_n = true;
		p.altrng_requested = false;  // No altitude filter
		in_region = false;
		break;
          case 'S':
		strcpy(p.satname, argv[i] + par_pos);
		len_satname = strlen(p.satname);
		p.satname_filter = true;
		break;
          case 't':
		p.delta_time = atof(argv[i] + par_pos);
		break;
          default:
		open_json();
		sprintf(errmsg, "Unrecognized command-line option '%s'", argv[i]);
  		close_stat_json(-2, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		exit(-2);
		break;
	}  // end switch
  }  // end if - for

  jd = JD_OFF + p.mjd;  // From MJD to JD

/* Use default TLE repository dir. ? */
  if ( p.use_deftledir )
	strcpy(tle_path, DEF_TLEDIR);

/* Figure out where the observer really is in Cartesian coordinates of date */
  earth_lat_alt_to_parallax(p.lat * DEG2RAD, p.ht_in_meters, &rho_cos_phi, &rho_sin_phi);
  observer_cartesian_coords(jd,
	p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, observer_loc);
  observer_cartesian_coords(jd + p.delta_time / SEC_IN_DAY,
	p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, observer_loc2);
 
  const short fname_len = 64;
  char *s, row[fname_len], norad_name[7], tle_list_file[200];
  bool input_list = false;
  FILE *lisfile, *ifile;

  open_json();


/* TLE file passed as first or last parameter? */
  int iarg = 1;
  if ( argv[1][0] == '-' )
	iarg = argc - 1;

  if ( *argv[iarg] == '@' ) {
	input_list = true;
	//tle_list_file = ++argv[1];
	sprintf(tle_list_file, "%s/%s", tle_path, ++argv[iarg]);
	if ( !(lisfile = fopen(tle_list_file, "rb")) ) {
		sprintf(errmsg, "Could not open input file %s", tle_list_file);
  		close_stat_json(2, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		exit(2);
	}
  }

/* HEALPix stuff */
  int order = 8, nside = 1 << order;
  long int id8nest;
  double theta;

/* A buffer to save satellite names and check for duplications in the various tle files */
  const int max_satbuff = 70000, extra_satbuff = 10000;
  int use_bufflen = 0, cur_bufflen = max_satbuff;
  double hareg;
  char *sbuff, *tle_file_name;
  bool single_sat = false, single_sat_found = false, read_tle_list = true;

  if ( (sbuff = (char *)malloc(max_satbuff)) == NULL ) {
	sprintf(errmsg, "Could not allocate required buffer memory");
  	close_stat_json(3, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
	exit(3);
  }

  sbuff[0] = '\0';

  p.tle_file_name = argv[iarg];

  if ( p.single_sat_n || p.single_sat_i ) {
	single_sat = true;
	p.geoloc_requested = true;
  }

/* Local Mean Sidereal Time & GMST */
  //double gmst;
  p.lmst = lmst_hr(jd, p.lon, &p.gmst);

/* For input geoloc or single satellite, use zenith coords if not given */
  if ( p.geoloc_reference || (single_sat && !in_region) ) {
	p.ra_deg = p.lmst * 15.;  
	p.de_deg = p.lat;
	target_ra = p.lmst / RAD2HRS;
	target_dec = p.lat * DEG2RAD;
	target_lon = p.lon * DEG2RAD;
  	target_lat = p.lat * DEG2RAD;
	hareg = 0.;
	add_geofieldsdesc_json();
  } else {
	target_ra = p.ra_deg * DEG2RAD;
	target_dec = p.de_deg * DEG2RAD;
	hareg = p.lmst - p.ra_deg / 15.;
  }

/* Alt, Az, Parang from Dec, HA, Lat */
  dechalat2alt(target_dec, hareg, p.lat, &p.alt, &p.az, &p.parang);

  add_params_json(p);

/* Sun position */
  Sun sun;

  lpsun_radec(jd, &sun.ra, &sun.dec);

/* Lon in [0, 360[ deg +W */
  sun.lon = fmod(p.gmst * 15. + 360. - sun.ra * RAD2DEG, 360.);
 
  sun.ha = p.lmst - sun.ra * RAD2HRS;
  dechalat2alt(sun.dec, sun.ha, p.lat, &sun.alt, &sun.az, &sun.parang);

  if ( p.geoloc_reference )
    sun.sep = skysep_h(target_lon, target_lat, sun.lon*DEG2RAD, sun.dec);
  else
    sun.sep = skysep_h(target_ra, target_dec, sun.ra, sun.dec);


/* Sun cartesian coordinates (in km) */ 
  sphrad2v(sun.ra, sun.dec, sun.d_km);

  for (i = 0; i < 3; i++)
	sun.d_km[i] *= AU_KM;


  sun.ra *= RAD2DEG;
  sun.dec *= RAD2DEG;

/* Lon in [-180, +180] deg +E */
  if ( sun.lon > 180. )
	sun.lon = 360. - sun.lon;
  else
	sun.lon *= -1;

#ifdef DEBUG
printf("JD %lf \n", jd);
printf("Region HA, LMST, AZ, Alt, PA: %lf %lf %lf %lf %lf (h, deg, deg, deg)\n", hareg, p.lmst, p.az, p.alt, p.parang);
printf("Sun RA, Dec, sep: %lf %lf  %lf  %lf (deg)\n", sun.ra, sun.dec, sun.lon, sun.sep);
printf("Sun HA, AZ, Alt, PA: %lf %lf %lf %lf (h, deg, deg, deg)\n", sun.ha, sun.az, sun.alt, sun.parang);
#endif
  add_sundata_json(sun);

/* Geodetic position */
  Geoloc geo;

/* Loop on list of tle files or simple open of single input file */
  while ( read_tle_list ) {

    char *lastchar;
    if ( !input_list ) {
	if ( strchr(tle_file_name, '/') == NULL )
	  sprintf(tle_path_file, "%s/%s", tle_path, argv[iarg]);
	else
	  strcpy(tle_path_file, tle_file_name);
	read_tle_list = false;
    } else {

	tle_file_name = row;
	if ( fgets(row, fname_len, lisfile) == NULL )
	  break;
	if ( row[0] == '#' )  /* Skip lines with comments */
	  continue;
	lastchar = tle_file_name + strlen(tle_file_name) - 1;  /* Drop linefeeds, CR, ... */
	if (*lastchar < 32)
	  *lastchar = (char) 0;
	if ( strlen(row) == 0 )  /* Skip empty lines */
	  continue;
	sprintf(tle_path_file, "%s/%s", tle_path, tle_file_name);
    }
 
    if ( !(ifile = fopen(tle_path_file, "rb")) ) {  /* Report error message and stop */
	sprintf(errmsg, "Could not open input file %s", tle_path_file);
  	close_stat_json(2, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
	exit(2);
    }


    if ( !fgets( line1, sizeof(line1), ifile) ) {
	sprintf(errmsg, "Could not read first TLE line from file %s", tle_path_file);
	close_stat_json(-2, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
	exit(-2);
    }

    strncpy(sat_name, line1, 23);
    trimend(sat_name);

    while ( fgets(line2, sizeof(line2), ifile) ) {

      tle_t tle;  /* Structure for two-line elements set for satellite */

      if ( !parse_elements(line1, line2, &tle) ) { /* TLE found */
/* First check if we have processed this sat already */
	strncpy(norad_name, line1+2, 5);  /* this is tle.norad_number */
	norad_name[5] = '\0';
//printf("\Satname:->%s\n", sat_name);
//printf("NORAD:->%s\n", norad_name);

	 char intl_desig[12] = {"none"};
	 if ( p.single_sat_n ) {
		int norad_name_n = strtol(norad_name, &endptr, 10);
		if ( norad_name_n == p.norad_n )
			single_sat_found = true;
	 } else if ( p.single_sat_i ) {
		if ( tle.intl_desig[0] != ' ' )
		  snprintf( intl_desig, 10, "%s%.2s-%s",
			(atoi( tle.intl_desig) > 57000 ? "19" : "20"),
			tle.intl_desig, tle.intl_desig + 2);

		if ( !memcmp(intl_desig, p.intl_desig, len_intl_desig) )
		  single_sat_found = true;
	} else if ( p.satname_filter ) {
		if ( memcmp(sat_name, p.satname, len_satname) )
		continue;
	}
 
/* Single satellite requested? */
	if ( single_sat && !single_sat_found )
		continue;
/* Add a , to make it unique, and skip any previously met satellite */
	strcat(norad_name, ",");
	if ( (s = strstr(sbuff, norad_name)) != NULL ) {
//printf("Found string '%s' at index = %ld ... skipping.\n", norad_name, s - sbuff);
		continue;
	} else {
		if ( use_bufflen >= cur_bufflen ) {
		  cur_bufflen += extra_satbuff;
//printf("\nRealloc extra bytes. Now %d \n", cur_bufflen);
		  if ( (sbuff = (char *)realloc(sbuff, cur_bufflen)) == NULL ) {
		    sprintf(errmsg, "Could not allocate required buffer memory");
		    close_stat_json(1, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		    exit(1);
		  }
		}
		strcat(sbuff, norad_name);
		use_bufflen += 6;
	}

	int is_deep = select_ephemeris(&tle),
	    is_in_sunlight;
//printf("is_deep: %d\n", is_deep);
	double sat_params[N_SAT_PARAMS], ang_sep, ang_sep1, d_ra, d_dec,
		ra, dec, ra1, dec1, sep_to_satellite, t_since,
		pos[3],  /* Satellite position vector */
		unused_delta2;

	t_since = (jd - tle.epoch) * 1440.;  /* From days to minutes */

	if ( is_deep ) {
		SDP4_init(sat_params, &tle);
		SDP4(t_since, &tle, sat_params, pos, NULL);
	} else {
		SGP4_init(sat_params, &tle);
		SGP4(t_since, &tle, sat_params, pos, NULL);
	}
	get_satellite_ra_dec_delta(observer_loc, pos, &ra, &dec, &sep_to_satellite);
/* For malformed or partial tle (e.g. Bepi-Colombo) this could happen */
	if ( isnan(ra) || isnan(dec) ) {
		sprintf(errmsg, "Could not decode TLE data");
		close_stat_json(2, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);
		exit(1);
	}

/* For single satellite also report its geo Lat, Lon, Alt, theta...
   but also computed when explicitly requested by the user (and added below)
*/
	if ( p.geoloc_requested )
		sat_geoLocation(p.gmst, pos,  &geo);

/* If altitude filter requested */
	if ( p.altrng_requested )
		if ( geo.alt < p.alt_min || geo.alt > p.alt_max )
			continue;

	epoch_of_date_to_j2000(jd, &ra, &dec);  /* Approx precession. Returned RA, Dec in radians. */

/* Compute position delta_time seconds later to
   1. check if object enters the requested region,
   2. compute speed/PA of motion (and then motion direction)

  Note: for geolocation only compute the reference position, not that at the second epoch
*/
	t_since += p.delta_time * 1440. / SEC_IN_DAY;

	if ( is_deep )
		SDP4(t_since, &tle, sat_params, pos, NULL);
	else
		SGP4(t_since, &tle, sat_params, pos, NULL);

	get_satellite_ra_dec_delta(observer_loc2, pos, &ra1, &dec1, &unused_delta2);
	epoch_of_date_to_j2000(jd, &ra1, &dec1);
/* Approx or precise separation? */
	if ( !p.haversine ) {
	  if ( p.geoloc_reference ) {
		d_ra = geo.lon * DEG2RAD - target_lon + TWOPI;
		if ( d_ra > PI )
		  d_ra -= TWOPI;
		d_dec = geo.lat - target_lat;
		ang_sep = sqrt(d_ra * d_ra + d_dec * d_dec) * RAD2DEG;
		ang_sep1 = 2 * ang_sep;  // dummy
	  } else {
		d_ra = (ra - target_ra + TWOPI);
		if ( d_ra > PI )
		  d_ra -= TWOPI;
		d_dec = dec - target_dec;
		ang_sep = sqrt(d_ra * d_ra + d_dec * d_dec) * RAD2DEG;
		d_ra = (ra1 - target_ra + TWOPI);
		if ( d_ra > PI )
		  d_ra -= TWOPI;
		d_dec = dec1 - target_dec;
		ang_sep1 = sqrt(d_ra * d_ra + d_dec * d_dec) * RAD2DEG;
	  }
	} else {
	  if ( p.geoloc_reference ) {
		ang_sep = skysep_h(target_lon, target_lat, geo.lon * DEG2RAD, geo.lat * DEG2RAD);
		ang_sep1 = 2 * ang_sep;  // dummy
//printf("\ntarget_lon, target_lat, geo.lon, geo.lat, ang_sep, search_radius: %lf %lf  %lf %lf %lf  %lf\n", target_lon, target_lat, geo.lon, geo.lat, ang_sep, p.search_radius);
	  } else {
		ang_sep = skysep_h(target_ra, target_dec, ra, dec);
		ang_sep1 = skysep_h(target_ra, target_dec, ra1, dec1);
	  }
	}

/* Check for single sat. or if in search area for both start and end epoch */
	if ( single_sat ||  /* for single satellite ignore region */
		p.satname_filter ||
		ang_sep <= p.search_radius || ang_sep1 <= p.search_radius ) {  /* good enough */
	  line1[16] = '\0';
	  double speed, posn_ang_of_motion, eray[2];

/* If satellite is in Sun light (this is preliminary. TODO to account for Sun angular size) */
	  is_in_sunlight = ! intersect_satsun_sphere(pos, sun.d_km, eray);
//printf("\n\n%s to Sun ray intersect Earth at %lf  %lf\n\n", sat_name, eray[0], eray[1]);
	  if ( p.in_sunlight_only && ! is_in_sunlight )
		continue;

	  if ( !p.single_sat_i ) {
		if ( tle.intl_desig[0] != ' ' )
		  snprintf( intl_desig, 10, "%s%.2s-%s",
			(atoi( tle.intl_desig) > 57000 ? "19" : "20"),
			tle.intl_desig, tle.intl_desig + 2);
	  }

	  d_ra = (ra1 - ra + TWOPI);
	  if ( d_ra > PI )
		d_ra -= TWOPI;
	  d_ra *= cos(dec);
	  d_dec = dec1 - dec;
	  posn_ang_of_motion = atan2(d_ra, d_dec);
	  if ( posn_ang_of_motion < 0. )
		posn_ang_of_motion += TWOPI;
//speed = sqrt(d_ra * d_ra + d_dec * d_dec) * RAD2DEG * 60 / p.delta_time;  // arcmin / s
	  speed = skysep_h(ra, dec, ra1, dec1) * 60 / p.delta_time;

/* HEALPix order 8 nested ID */
	  theta = (PI/2. - dec);
	  ang2pix_nest(nside, theta, ra, &id8nest);

	  if ( !p.info_only ) {
		if ( n_sats_found > 0 )
		  printf(", ");
		else {
		  add_fieldsdesc_json();
		  printf(" \"satellites\": [");
		}

		printf("{\"name\": \"%s\", \"intl_desig\": \"%s\", \"norad_n\": %d, ",
			sat_name, intl_desig, tle.norad_number);

		if ( p.geoloc_requested )
			add_satlatlon_json(geo);

		printf("\"data\": [%8.4lf,%8.4lf,%8.4lf,%8.4lf,%9.2lf,%6.4lf,%3d,%6.3lf,%d,%ld]}",
			ra * RAD2DEG, dec * RAD2DEG, ra1 * RAD2DEG, dec1 * RAD2DEG,
			sep_to_satellite, ang_sep,
			(int)(posn_ang_of_motion * RAD2DEG), speed, is_in_sunlight, id8nest);
/* Speed is displayed in arcminutes/second (== degrees/minute) */
	  }
	  n_sats_found++;

	  if ( is_in_sunlight )
		n_sats_in_sunlight++;

	  if ( single_sat_found || n_sats_found == p.max_sats ) {
		fclose(ifile);
		goto end_checking;
	  }
	} // end if ang_sep < p.search_radius
      }  // end TLE found

      strcpy(line1, line2);
      if ( line2[0] != '1' && line2[0] != '2' ) {
	  strncpy(sat_name, line2, 23);
	  trimend(sat_name);
      }
    }  // end while fgets(line2 ...

    fclose(ifile);
  }  // end while read_tle_list


end_checking:
  if ( single_sat_found )
	n_sats = 1;
  else
	n_sats = n_sats_found;  /* This is TODO using HEALPix IDs */

  free(sbuff);

  if ( !p.info_only && n_sats_found > 0 )
	printf("],");

/* For single satellite request report failure but do not rise an error flag */
  if ( single_sat && !single_sat_found )
	sprintf(errmsg, "Requested satellite not found in TLE file");

  close_stat_json(status, errmsg, n_sats_found, n_sats, n_sats_in_sunlight);

  return(0);
}
