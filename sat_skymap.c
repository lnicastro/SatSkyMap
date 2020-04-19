/*
   sat_skymap.c

  Purpose:
    Find which satellites are within a given angular separation (approx)
    of a given RA/Dec (J2000), as seen from a given point.
    Single satelletite info, regardless of its position can be requested, too.

  Procedure:
    The code reads in a TLE file (name provided as the first command-line
    argument) with details of the observer position, search radius, date/time,
    and RA/Dec provided on the command line.

  Input Parameters:
    -d CalendarDate	Calendar date (UTC) of interest (in the form yyyy-mm-ddThh:mm:ss[.sss])
    -i sat_intnlname	Single satellite selection via its international designator (region ignored)
    -j MJD		Modified Julian Date of interest (ignored if Calendar Date given)
    -l Lat,Lon,Height	Observation site (comma separated data with no spaces)
    -n MaxSats		Maximum number of satellites to return (def. 1000)
    -p RA,Dec		J2000 sky coordinates of region to check
    -r radius		Region radius centered at the given coords
    -s sat_norad_n	Single satellite selection via its NORAD number (region ignored)
    -t deltat		Second epoch delta time (s, def. 1)
    -D DirTLEs		Directory with the repository of TLE files
  Switches:
    -h			print this help
    -H			Compute sky separation via Haversine formula rather than cartesian triangle (suggested! Default?)
    -I			Information about the returned data and number of satellites found
  			('satellites' object not returned)

  Output:
    A json string with information about the satellites in the FoV (see below).

  Example (La Silla coords):

    ./sat_skymap default.tle -l-29.25627,-70.73805,2400 -j58861.5 -p90.5,-30.3 -r20
    ./sat_skymap default.tle -l-29.25627,-70.73805,2400 -d2020-01-13T12:00:00 -p90.5,-30.3 -r20

  Scan the TLE element file 'default.tle' for satellites visible from
  (lat, lon, height) = (-29.25627, -70.7380, 2400),
  on MJD 58861.5 (UTC: 2020-01-13 12h = 2020-01-13T12:00:00),
  at (RA, Dec) = 90.5, +30.3 (deg), within a 20-deg search radius.
  Positions at given JD + deltat (def. 1) s is also computed and returned.
  Additional computed info:
    Local Mean Sidereal Time
    region Az, Alt, Parang
    Sun RA/Dec, Alt/Az coordinates.

  The output looks like this:

  {
  "swinfo": {"name": "sat_skymap", "author": "L. Nicastro @ INAF-OAS", "date": "2020-02-26", "version": "0.1b"},
  "input_params": {"tle_file": "default.tle", "location": ["lat":-29.2563, "long": -70.7381, "height":  2400.0],
  "region": {"ra":  90.5000, "dec":-30.3000, "radius": 20.0000, "lmst": 14.7803, "az": 222.1310, "alt":-14.4561, "parang": 137.324},
  "notes": "All coordinates and radius in degrees. LMST in hrs.",
  "mjd": 58861.50000, "epoch_UTC": "2020-01-13T12:00:00", "delta_time_s": 1, "max_sats": 1000},
  "sun": {"ra":294.566, "dec":-21.517, "az": 101.818, "alt": 24.735, "parang":-113.376, "separation_deg":123.255},
  "data_fields": {"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "HPXID_8"],
    "desc": ["RA Tinit", "Dec Tinit", "RA Tend", "Dec Tend", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "HEALPix order 8 nested schema ID"],
    "type": ["double", "double", "double", "double", "double", "float", "float", "float", "int"],
    "unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", ""]},
  "satellites": [{"name": "STARLINK-23", "intl_desig": "1919-029C ", "norad_n": 44237,
    "data": [ 77.5770, -17.4894,  77.6262, -17.4900,   7660.9,18.1955,  90,  2.818,339078]},


  { ... }}],
  "status": 0, "errmsg": "", "n_sats_found": 6, "n_sats": 6
  }

  Where:
    Distance: distance to satellite in km,
    Separation: angular separation in degrees from the search point,
    PA: position angle of motion,
    Speed: apparent angular rate of motion in arcminutes/second (or degrees/minute),
    HPXID_8: HEALPix order 8 nested schema ID

  Example (ISS position):
    ./sat_skymap stations.txt -H -D/usr/local/TLErepo -l44.52804,11.33715,23.5 -j58941.71931 -s 25544

  Output:

  {
  "swinfo": {"name": "sat_skymap", "author": "L. Nicastro @ INAF-OAS", "date": "2020-04-06", "version": "0.1d"},
  "input_params": {"tle_file": "stations.txt", "location": {"lat": 44.5280, "long":  11.3371, "height":    23.5},
  "region": {"ra": 101.7991, "dec": 44.5280, "radius": 20.0000, "lmst":  6.7866, "az": 190.0614, "alt": 14.4649, "parang":   8.294},
  "notes": "All coordinates and radius in degrees. LMST in hrs.",
  "mjd": 58941.71931, "epoch_UTC": "2020-04-02T17:15:48", "delta_time_s": 1, "max_sats": 1000},
  "sun": {"ra": 12.353, "dec":  5.298, "az": 273.396, "alt":  4.107, "parang":  45.619, "separation_deg": 82.530},
  "data_fields": {"name": ["RA_start", "Dec_start", "RA_end", "Dec_end", "Distance", "Separation", "PA", "Speed", "HPXID_8"],
    "desc": ["RA T_ini", "Dec T_ini", "RA T_end", "Dec T_end", "distance to sat.", "angular separation", "position angle", "apparent angular rate of motion", "HEALPix order 8 nested schema ID"],
    "type": ["double", "double", "double", "double", "double", "float", "float", "float", "int"],
    "unit": ["deg", "deg", "deg", "deg", "km", "deg", "deg", "arcmin/s", ""]},
  "satellites": [{"name": "ISS (ZARYA)", "intl_desig": "1998-067A ", "norad_n": 25544,
    "data": [249.7761,-45.8183,249.8236,-45.8294, 12257.7,101.5945,108, 2.096,675636]}], "status": 0, "errmsg": "", "n_sats_found": 1, "n_sats": 1}



   LN @ INAF-OAS, Jan 2020.  Last change: 16/04/2020
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

  typedef struct myParams {
	char *tle_file_name;
	double lat;
	double lon;
	double ht_in_meters;
	double ra_deg;
	double de_deg;
	double search_radius;
	double az;
	double alt;
	double parang;
	double lmst;
	double mjd;
	char date[24];
	char intl_desig[12];
	int delta_time;
	int max_sats;
	int norad_n;
	bool haversine;
	bool info_only;
	bool single_sat_i;
	bool single_sat_n;
  } Params;

  typedef struct mySun {
	double ra;
	double dec;
	double az;
	double alt;
	double parang;
	double sep;
  } Sun;

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
  "  -d CalendarDate	Calendar date (UTC) of interest (in the form yyyy-mm-ddThh:mm:ss[.sss])\n"
  "  -i sat_intnlname	Single satellite selection via its international designator (region ignored)\n"
  "  -j MJD		Modified Julian Date of interest (ignored if Calendar Date given)\n"
  "  -l lat,lon,height	Observation site (comma separated data with no spaces)\n"
  "  -n MaxSats		Maximum number of satellites to return (def. 1000)\n"
  "  -p Ra,Dec		J2000 sky coordinates of region to check\n"
  "  -r radius		Region radius centered at the given coords\n"
  "  -s sat_norad_n	Single satellite selection via its NORAD number (region ignored)\n"
  "  -t deltat		Second epoch delta time (s, def. 1)\n"
  "  -D DirTLEs		Directory with the repository of TLE files\n\n"
  "Switches:\n"
  "  -h			print this help\n"
  "  -H			Compute sky separation via Haversine formula rather than cartesian triangle (suggested! Default?)\n"
  "  -I			Information about the returned data and number of satellites found\n"
  " 			('satellites' object not returned)\n\n"

  "Example usage:\n"
  "  %s default.tle -l-29.25627,-70.73805,2400 -p90.5,-30.3 -j58861.5 -r20\n"
  "  %s stations.txt -H -D/usr/local/TLErepo -l44.52804,11.33715,23.5 -j58941.71931 -s 25544\n", progname, progname, progname);

// exit(0);
}


void open_json() {
  printf("{"
    "\"swinfo\": {"
    "\"name\": \"%s\", \"author\": \"%s\", \"date\": \"%s\", \"version\": \"%s\"}, ",
    progname, progauthor, progdate, progversion);
}

void add_params_json(Params p) {
  printf("\"input_params\": {"
    "\"tle_file\": \"%s\", \"location\": {\"lat\":%8.4lf, \"long\":%9.4lf, \"height\":%8.1lf}, \"region\": {\"ra\":%9.4lf, \"dec\":%8.4lf, \"radius\":%8.4lf, \"lmst\":%8.4lf, \"az\":%9.4lf, \"alt\":%8.4lf, \"parang\":%8.3lf}, \"notes\": \"All coordinates and radius in degrees. LMST in hrs.\", "
    "\"mjd\": %11.5lf, \"epoch_UTC\": \"%s\", \"delta_time_s\": %d, \"max_sats\": %d}, ",
    p.tle_file_name, p.lat, p.lon, p.ht_in_meters, p.ra_deg, p.de_deg, p.search_radius, p.lmst, p.az, p.alt, p.parang, p.mjd, p.date, p.delta_time, p.max_sats);
}

void add_sundata_json(Sun sun) {
  printf("\"sun\": {"
    "\"ra\":%7.3lf, \"dec\":%7.3lf, \"az\":%8.3lf, \"alt\":%7.3lf, \"parang\":%8.3lf, \"separation_deg\":%7.3lf}, ",
    sun.ra, sun.dec, sun.az, sun.alt, sun.parang, sun.sep);
}

void add_fieldsdesc_json() {
  printf(
    "\"data_fields\": {\"name\": [\"RA_start\", \"Dec_start\", \"RA_end\", \"Dec_end\", \"Distance\", \"Separation\", \"PA\", \"Speed\", \"HPXID_8\"], "
    "\"desc\": [\"RA T_ini\", \"Dec T_ini\", \"RA T_end\", \"Dec T_end\", \"distance to sat.\", \"angular separation\", \"position angle\", \"apparent angular rate of motion\", \"HEALPix order 8 nested schema ID\"], "
    "\"type\": [\"double\", \"double\", \"double\", \"double\", \"double\", \"float\", \"float\", \"float\", \"int\"], "
    "\"unit\": [\"deg\", \"deg\", \"deg\", \"deg\", \"km\", \"deg\", \"deg\", \"arcmin/s\", \"\"]},");
}

void add_satlatlon_json(double lat, double lon, double alt) {
  printf( "\"geoloc\": {\"lat\":%7.3lf, \"lon\":%7.3lf, \"alt\":%9.2lf},", lat, lon, alt);
}

void close_stat_json(int status, char *errmsg, int n_sats_found, int n_sats) {
  printf(" \"status\": %d, \"errmsg\": \"%s\", \"n_sats_found\": %d, \"n_sats\": %d"
	"}\n", status, errmsg, n_sats_found, n_sats);
}


/* GMST in hours. Reference:  The 1992 Astronomical Almanac, page B6. */
/*
static inline double myThetaG( double jd)
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


/* Compute satellite geodetic Lon, Lat and Alt position given its ECI position and GMST.
 * Returned values are deg, deg, km.
 * */
void sat_geoLocation(double gmst, double *pos,  double *lon, double *lat, double *alt) {
  double r, e2, phi, sinphi, c;
  const double xkmper = 6.378137e3;      /* WGS 84 Earth radius km */
  const double f = 3.35281066474748e-3;  /* Flattening factor */

  double theta = atan2(pos[1], pos[0]) * RAD2DEG;  /* degrees */
  //lon = (theta - gmst) % TWOPI;
  *lon = (theta - gmst * 15.);
//if (*lon > 360.)
  if (*lon > 180)
	*lon -= 360.;
//else if (*lon < 360.)
  else if (*lon < -180.)
	*lon += 360.;

  r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
  e2 = f * (2. - f);
  *lat = atan2(pos[2], r);

  do {
	phi = *lat;
	sinphi = sin(phi);
	c = 1. / sqrt(1. - e2 * sinphi*sinphi);
	*lat = atan2(pos[2] + xkmper * c * e2 * sinphi, r);
  } while (fabs(*lat - phi) >= 1e-10);

  *alt = r/cos(*lat) - xkmper * c;  /* kilometers */

  *lat *= RAD2DEG;  /* degrees */
   if (*lat > 90.)
	*lat -= 360.;
}


int main(int argc, char **argv)
{
  char tle_path[128] = {"."}, errmsg[128] = {""}, line1[100], line2[100], sat_name[25], tle_path_file[200],
	opt, *endptr;
  double jd, target_ra, target_dec, rho_sin_phi, rho_cos_phi, observer_loc[3], observer_loc2[3];
  int i, par_pos, header_line_shown = 0;
  int status = 0, n_sats_found = 0, n_sats = 0;
  bool in_region = false;  // used for single satellite request


  if ( argc < 2 ) {
	open_json();
	sprintf(errmsg, "Example usage:\n %s default.tle -l-29.25627,-70.73805,2400 -p90.5,-30.3 -j58861.5 -r20", progname);
 	close_stat_json(1, errmsg, n_sats_found, n_sats);
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
  p.intl_desig[0] = '\0';
  p.norad_n = 0;


// Note: no check on validity of input paramters!
//
  for( i = 1; i < argc; i++ )
    if( argv[i][0] == '-') {
//printf("\n In %d '%s'  %s***\n", i, &(argv[i][2]), argv[i]);
		opt = argv[i][1];
		par_pos = 2;
		if ( argv[i][2] == 0 ) {  // blank after option
			//if ( argv[i+1][0] != '-' )  // No, because of possible negative numbers after the param
		  ++i;
		  par_pos = 0;
		}
//printf("Read %d  %d %s  %s***\n", i, par_pos, &opt, argv[i]);
	switch( opt )
	{
	  /* Switches */
 	  case 'h':
		Usage();
		exit(0);
 	  case 'H':
		p.haversine = true;
		i--;
		break;
 	  case 'I':
		p.info_only = true;
		i--;
		break;
	  /* Parameters */
          case 'd':
		strcpy(p.date, argv[i] + par_pos);
		date2mjd(argv[i] + par_pos, &p.mjd);
		break;
 	  case 'D':
		strcpy(tle_path, argv[i] + par_pos);
		break;
          case 'i':
		strcpy(p.intl_desig, argv[i] + par_pos);
		p.single_sat_i = true;
		break;
          case 'j':
		p.mjd = atof(argv[i] + par_pos);
		strcpy(p.date, mjd2date(p.mjd));
		break;
          case 'l':
		sscanf(argv[i] + par_pos, "%lf,%lf,%lf", &p.lat, &p.lon, &p.ht_in_meters);
		break;
          case 'n':
		sscanf(argv[i] + par_pos, "%d", &p.max_sats);
		break;
          case 'p':
		sscanf(argv[i] + par_pos, "%lf,%lf", &p.ra_deg, &p.de_deg);
		if ( p.ra_deg < 0. || p.ra_deg > 360. )  {
		  open_json();
		  sprintf(errmsg, "RA must be in the range 0, 360 degrees. Read '%lf'", p.ra_deg);
  		  close_stat_json(-1, errmsg, n_sats_found, n_sats);
		  exit(-1);
		}
		if ( p.de_deg < -90. || p.de_deg > 90. )  {
		  open_json();
		  sprintf(errmsg, "Dec must be in the range -90, +90 degrees. Read '%lf'", p.de_deg);
  		  close_stat_json(-1, errmsg, n_sats_found, n_sats);
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
		break;
          case 't':
		p.delta_time = atof(argv[i] + par_pos);
		break;
          default:
		open_json();
		sprintf(errmsg, "Unrecognized command-line option '%s'", argv[i]);
  		close_stat_json(-2, errmsg, n_sats_found, n_sats);
		exit(-2);
		break;
	}  // end switch
  }  // end if - for

  jd = JD_OFF + p.mjd;  // From MJD to JD

/* Figure out where the observer _really_ is in Cartesian coordinates of date */
  earth_lat_alt_to_parallax(p.lat * DEG2RAD, p.ht_in_meters, &rho_cos_phi, &rho_sin_phi);
  observer_cartesian_coords(jd,
	p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, observer_loc);
  observer_cartesian_coords(jd + TIME_EPSILON * p.delta_time,
	p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, observer_loc2);
 
  const short fname_len = 64;
  char *s, row[fname_len], norad_name[6], tle_list_file[200];  // const char *tle_list_file;
  bool input_list = false;
  FILE *lisfile, *ifile;

  open_json();

  if ( *argv[1] == '@' ) {
	input_list = true;
	//tle_list_file = ++argv[1];
	sprintf(tle_list_file, "%s/%s", tle_path, ++argv[1]);
	if ( ! (lisfile = fopen(tle_list_file, "rb")) ) {
		sprintf(errmsg, "Couldn't open input file %s", tle_list_file);
  		close_stat_json(2, errmsg, n_sats_found, n_sats);
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
  double hareg, hasun;
  char *sbuff, *tle_file_name;
  bool single_sat_found = false, read_tle_list = true;

  if ( (sbuff = (char *)malloc(max_satbuff)) == NULL ) {
	sprintf(errmsg, "Could not allocate required buffer memory");
  	close_stat_json(3, errmsg, n_sats_found, n_sats);
	exit(3);
  }

  sbuff[0] = '\0';

  p.tle_file_name = argv[1];

  double gmst;
  p.lmst = lst_hr(jd, p.lon, &gmst);  // Local Mean Sidereal Time
//printf("lst_hr gmst: %lf\n", gmst);
//double tmp = myThetaG(jd);
//printf("ThetaG: %lf\n", tmp);

// For single satellite, use zenith coords if not given
  if ( (p.single_sat_n || p.single_sat_i) && !in_region ) {
	p.ra_deg = p.lmst * 15.;  
	p.de_deg = p.lat;
	hareg = 0.;
	//p.az = 0;
	//p.alt = 90.;
  } else
	hareg = p.lmst - p.ra_deg / 15.;

  target_ra = p.ra_deg * DEG2RAD;
  target_dec = p.de_deg * DEG2RAD;

  dechalat2alt(target_dec, hareg, p.lat, &p.alt, &p.az, &p.parang);

  add_params_json(p);

/* Sun position */

  Sun sun;

  lpsun_radec(jd, &sun.ra, &sun.dec);
  sun.sep = skysep_h(target_ra, target_dec, sun.ra, sun.dec);
 
  hasun = p.lmst - sun.ra/HRS2RAD;
  dechalat2alt(sun.dec, hasun, p.lat, &sun.alt, &sun.az, &sun.parang);

  sun.ra *= RAD2DEG;
  sun.dec *= RAD2DEG;

#ifdef DEBUG
printf("JD %lf \n", jd);
printf("Region HA, LMST, AZ, Alt, PA: %lf %lf %lf %lf %lf (h, deg, deg, deg)\n", hareg, lmst, p.az, p.alt, p.parang);
printf("Sun RA, Dec, sep: %lf %lf  %lf (deg)\n", sun.ra, sun.dec, sun.sep);
printf("Sun HA, AZ, Alt, PA: %lf %lf %lf %lf (h, deg, deg, deg)\n", hasun, sun.az, sun.alt, sun.parang);
#endif

  add_sundata_json(sun);


/* Loop on list of tle files or simple open of single input file */
  while ( read_tle_list ) {

	char *lastchar;
	if ( !input_list ) {
		if ( strchr(tle_file_name, '/') == NULL )
		  //strcat(tle_path, "/");
		  sprintf(tle_path_file, "%s/%s", tle_path, argv[1]);
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
		 //trimend(row);
		if (*lastchar < 32)
		  *lastchar = (char) 0;
		if ( strlen(row) == 0 )  /* Skip empty lines */
		  continue;
		sprintf(tle_path_file, "%s/%s", tle_path, tle_file_name);
//printf("tle_path_file=%s\n", tle_path_file);
	}
 
	if ( !(ifile = fopen(tle_path_file, "rb")) ) {  /* Report error message but do not stop */
		sprintf(errmsg, "Couldn't open input file %s", tle_path_file);
  		close_stat_json(0, errmsg, n_sats_found, n_sats);
		continue;
	}


 	if ( !fgets( line1, sizeof(line1), ifile) ) {
		sprintf(errmsg, "Couldn't read first TLE line from file %s", tle_path_file);
		close_stat_json(-2, errmsg, n_sats_found, n_sats);
		exit(-2);
	}

//line1[24] = '\0';
	strncpy(sat_name, line1, 23);
	trimend(sat_name);

    while ( fgets(line2, sizeof(line2), ifile) ) {

      tle_t tle;     /* Structure for two-line elements set for satellite */

      if ( !parse_elements(line1, line2, &tle) )    /* got a TLE */
      {
/* First check if we have processed this sat already */
	 strncpy(norad_name, line1+2, 5);  /* this is tle.norad_number */
	 norad_name[5] = '\0';
//printf("NORAD:->%s\n", norad_name);

	  char intl_desig[12] = {"none"};
	  if ( p.single_sat_n ) {
		int norad_name_n = strtol(norad_name, &endptr, 10);
		if ( norad_name_n == p.norad_n )
			single_sat_found = true;
	  } else if ( p.single_sat_i ) {

		if ( tle.intl_desig[0] != ' ' ) {
			short launch_year;
			strncpy(intl_desig, tle.intl_desig, 2);
			sscanf(intl_desig, "%hd", &launch_year);
			if ( launch_year > 99 )
			  strcpy(intl_desig, "20");
			else
			  strcpy(intl_desig, "19");
			strncat(intl_desig, tle.intl_desig, 2);
			strcat(intl_desig, "-");
			strcat(intl_desig, &line1[11]);
		}

		if ( intl_desig == p.intl_desig )
			single_sat_found = true;
	    }
 
// Single satellite requested?
	 if ( (p.single_sat_i || p.single_sat_n) && !single_sat_found )
		continue;

	 if ( (s = strstr(sbuff, norad_name)) != NULL ) {
//printf("Found string at index = %d ... skipping.\n", s - sbuff);
		continue;
	 } else {
		if ( use_bufflen >= cur_bufflen ) {
		  cur_bufflen += extra_satbuff;
//printf("\nRealloc extra bytes. Now %d \n", cur_bufflen);
		  if ( (sbuff = (char *)realloc(sbuff, cur_bufflen)) == NULL ) {
		    sprintf(errmsg, "Could not allocate required buffer memory");
		    close_stat_json(1, errmsg, n_sats_found, n_sats);
		    exit(1);
		  }
		}
		strcat(sbuff, norad_name);
		use_bufflen += 5;
	 }

	int is_deep = select_ephemeris(&tle);
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
// For malformed or partial tle (eg Bepi-Colombo) this could happen
	if ( isnan(ra) || isnan(dec) ) {
		sprintf(errmsg, "Could not decode properly TLE data");
		close_stat_json(2, errmsg, n_sats_found, n_sats);
		exit(1);
	}

/* For single satellite also report its Lat, Lon */
	if ( p.single_sat_i || p.single_sat_n ) {
		double sat_lat, sat_lon, sat_alt;
		sat_geoLocation(gmst, pos,  &sat_lon, &sat_lat, &sat_alt);
//printf("\n sat_lat, sat_lon, sat_alt: %lf, %lf, %lf, \n", sat_lat, sat_lon, sat_alt);
		add_satlatlon_json(sat_lat, sat_lon, sat_alt);

	}

	epoch_of_date_to_j2000(jd, &ra, &dec);  /* Approx precession. Returned RA, Dec in radians. */

/* Compute position delta_time seconds later to
 * 1. check if object enters the requested region,
 * 2. compute speed/PA of motion (and then motion direction)
 */
	t_since += TIME_EPSILON * p.delta_time * 1440.;

	if ( is_deep )
             SDP4(t_since, &tle, sat_params, pos, NULL);
	else
             SGP4(t_since, &tle, sat_params, pos, NULL);

	get_satellite_ra_dec_delta(observer_loc2, pos, &ra1, &dec1, &unused_delta2);
	epoch_of_date_to_j2000(jd, &ra1, &dec1);

/* Approx or precise separation? */
	if ( ! p.haversine ) {
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
	} else {
	   ang_sep = skysep_h(target_ra, target_dec, ra, dec);
	   ang_sep1 = skysep_h(target_ra, target_dec, ra1, dec1);
	}


/* Check for single sat. or if in search area for both start and end epoch */
	if ( p.single_sat_i || p.single_sat_n ||  /* for single satellite ignore region */
	     ang_sep < p.search_radius || ang_sep1 < p.search_radius )  /* good enough */
	{
            line1[16] = '\0';
            double speed, posn_ang_of_motion;

	    if ( ! p.single_sat_i ) {
	    //char intl_desig[12] = {"none"};

	      if ( tle.intl_desig[0] != ' ' ) {
	        short launch_year;
	        strncpy(intl_desig, tle.intl_desig, 2);
	        sscanf(intl_desig, "%hd", &launch_year);
	        if ( launch_year > 99 )
		  strcpy(intl_desig, "20");
	        else
		  strcpy(intl_desig, "19");
	        strncat(intl_desig, tle.intl_desig, 2);
	        strcat(intl_desig, "-");
	        strcat(intl_desig, &line1[11]);
	      }
	    }

	    if ( ! header_line_shown ) {
		if ( ! p.info_only )
		  printf(" \"satellites\": [");

		header_line_shown = 1;
	    }

            d_ra = (ra1 - ra + TWOPI);
            if ( d_ra > PI )
               d_ra -= TWOPI;
            d_ra *= cos(dec);
            d_dec = dec1 - dec;
            posn_ang_of_motion = atan2(d_ra, d_dec);
            if( posn_ang_of_motion < 0. )
               posn_ang_of_motion += TWOPI;
//speed = sqrt(d_ra * d_ra + d_dec * d_dec) * RAD2DEG * 60 / p.delta_time;  // arcmin / s
	    speed = skysep_h(ra, dec, ra1, dec1) * 60 / p.delta_time;

/* HEALPix order 8 nested ID */
	    theta = (PI/2. - dec);
	    ang2pix_nest(nside, theta, ra, &id8nest);
//printf( "RA, Dec, ID: %lf %lf  %ld\n", ra, dec, id8nest);

	    if ( !p.info_only ) {
		if ( n_sats_found > 0 )
		  printf(", ");
		else
		  add_fieldsdesc_json();

//if ( single_sat_found || (!p.single_sat_n && !p.single_sat_i) )
		printf( "{\"name\": \"%s\", \"intl_desig\": \"%s\", \"norad_n\": %d, \"data\": [%8.4lf,%8.4lf,%8.4lf,%8.4lf,%9.2lf,%6.4lf,%3d,%6.3lf,%ld]}",
			sat_name, intl_desig, tle.norad_number, ra * RAD2DEG, dec * RAD2DEG, ra1 * RAD2DEG, dec1 * RAD2DEG,
			sep_to_satellite, ang_sep,
			(int)(posn_ang_of_motion * RAD2DEG), speed, id8nest);
/* Speed is displayed in arcminutes/second, or in degrees/minute */
	    }
	    n_sats_found++;

	    if ( single_sat_found || n_sats_found == p.max_sats ) {
		fclose(ifile);
		goto end_checking;
	    }
          }  // end if ang_sep < p.search_radius
        }  // end got a TLE
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
	n_sats = n_sats_found;  // This is TODO using HEALPix IDs

  free(sbuff);

  if ( ! p.info_only && n_sats_found > 0 )
	printf("],");

/* For single satellite request report failure but do not rise an error flag */
  if ( (p.single_sat_n || p.single_sat_i) && ! single_sat_found )
	 sprintf(errmsg, "Requested satellite not found in TLE file");

  close_stat_json(status, errmsg, n_sats_found, n_sats);

  return(0);
}
