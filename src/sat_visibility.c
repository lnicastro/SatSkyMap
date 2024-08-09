/*
   sat_visibility.c

  Purpose:
    Find interval(s) of a satellite visibility (sunlit) at a given site given its NORAD number or international designator.

  Procedure:
    The code reads in a TLE file (name provided as the first or last command-line
    argument) with details of the observer position and date/time interval in UT.

  Input Parameters:
    -d CalendarDates	Calendar date (UTC) range of interest (in the form yyyy-mm-ddThh:mm:ss[.sss],yyyy-...)
    -D DirTLEs		Directory with the repository of TLE files (def. ./; ignore -T option)
    -i SatIntnlName	Satellite selection via its international designator
    -j MJDs		Modified Julian Date range of interest (MJDstart,MJDend - ignored if Calendar Date given)
    -l Lat,Lon,Alt	Geodetic observing site (comma separated data with no spaces)
    -r radius		Region search radius - centred at the site zenith (degrees; def. 90)
    -s SatNorad_n	Satellite selection via its NORAD number
    -t deltaT		Visibility check delta time (seconds; def. 3)

  Switches:
    -h			Print this help
    -N			Skip visibility check until Sun is set - if it is up at start date
    -T			Use default repository directory for TLE files (def. ./; see sat_skymap_def.h) 

  Output:
    A json string with information about the satellites visibility interval (if any).


  Example 1 (La Silla coords, MJD input, ISS position, check every 10 s):
    ./sat_visibility stations.txt -T -t 10 -l-29.25627,-70.73805,2400 -j60436.54,60436.55 -s 25544

  Example 2 (Bologna, Calendar dates, ISS position):
    ./sat_visibility stations.txt -T -l44.52324,11.33914,23 -d2024-05-11T20:00:00,2024-05-11T22:00:00 -s 25544


  For the latter case, the output will look like this:


  "swinfo": {
    "name": "sat_visibility",
    "author": "L. Nicastro @ INAF-OAS",
    "date": "2024-05-13",
    "version": "0.1a"
  },
  "geoloc_fields": {
    "lat": {
      "desc": "Geodetic Latitude",
      "unit": "deg"
    },
    "lon": {
      "desc": "Geodetic Longitude",
      "unit": "deg"
    },
    "alt": {
      "desc": "Geodetic Altitude",
      "unit": "km"
    },
    "theta": {
      "desc": "Equatorial angle (Lon + GMST = RA)",
      "unit": "deg"
    }
  },
  "input_params": {
    "tle_file": "stations.txt",
    "location": {
      "lat": 44.5232,
      "lon": 11.3391,
      "alt": 23
    },
    "region": {
      "ra": 0,
      "dec": 44.5232,
      "radius": 90,
      "lmst": 0,
      "az": 0,
      "alt": 90,
      "parang": 180
    },
    "mjd_range": [
      60441.83333,
      60441.91667
    ],
    "date_utc_range": [
      "2024-05-11T20:00:00",
      "2024-05-11T22:00:00"
    ],
    "gmst": 0,
    "delta_time_s": 3,
    "notes": "All coordinates and radius in degrees. GMST, LMST in hrs."
  },
  "sat_fields": {
    "desc": [
      "Date start/end sunlit",
      "MJD start/end sunlit",
      "RA J2000",
      "Dec J2000",
      "Azimuth",
      "Altitude",
      "Distance to sat.",
      "Angular separation",
      "Position angle",
      "Apparent angular rate of motion"
    ],
    "type": [
      "datetime",
      "double",
      "double",
      "double",
      "double",
      "double",
      "double",
      "float",
      "float",
      "float"
    ],
    "unit": [
      "",
      "MJD",
      "deg",
      "deg",
      "deg",
      "deg",
      "km",
      "deg",
      "deg",
      "arcmin/s"
    ]
  },
  "satellite": {
    "name": "ISS (ZARYA)",
    "intl_desig": "1998-067A",
    "norad_n": 25544,
    "geoloc_start": {
      "lat": 35.325,
      "lon": -12.141,
      "alt": 416.35,
      "theta": 169.37
    },
    "sun_start": {
      "ra": 49.218,
      "dec": 18.171,
      "az": 323.407,
      "alt": -19.074,
      "lon": -132.293,
      "parang": 26.572,
      "separation_deg": 109.074
    },
    "sat_start": {
      "Date": "2024-05-11T20:45:33",
      "MJD": 60441.86496,
      "ra": 116.1688,
      "dec": -13.3338,
      "az": 250.9996,
      "alt": 0.0387,
      "Distance": 2338.45,
      "Separation": 89.9613,
      "PA": 25,
      "Speed": 4.072
    },
    "geoloc_end": {
      "lat": 51.698,
      "lon": 36.589,
      "alt": 418.39,
      "theta": -139.335
    },
    "sun_end": {
      "ra": 49.225,
      "dec": 18.173,
      "az": 325.732,
      "alt": -20.13,
      "lon": -134.852,
      "parang": 24.994,
      "separation_deg": 110.13
    },
    "sat_end": {
      "Date": "2024-05-11T20:55:47",
      "MJD": 60441.87207,
      "ra": 307.4248,
      "dec": 23.7787,
      "az": 57.8708,
      "alt": 2.0926,
      "Distance": 2127.44,
      "Separation": 87.9074,
      "PA": 25,
      "Speed": 4.072
    },
    "geoloc_mid": {
      "lat": 46.459,
      "lon": 8.568,
      "alt": 417.95,
      "theta": -168.634
    },
    "sun_mid": {
      "ra": 49.221,
      "dec": 18.172,
      "az": 324.567,
      "alt": -19.612,
      "lon": -133.577,
      "parang": 25.788,
      "separation_deg": 109.612
    },
    "sat_mid": {
      "Date": "2024-05-11T20:50:41",
      "MJD": 60441.87207,
      "ra": 133.2947,
      "dec": 60.2011,
      "az": 315.7907,
      "alt": 51.5948,
      "Distance": 523.38,
      "Separation": 38.4052,
      "PA": 25,
      "Speed": 4.072
    },
    "geoloc_max": {
      "lat": 46.851,
      "lon": 9.696,
      "alt": 418,
      "theta": -167.447
    },
    "sun_max": {
      "ra": 49.222,
      "dec": 18.172,
      "az": 324.62,
      "alt": -19.636,
      "lon": -133.635,
      "parang": 25.752,
      "separation_deg": 109.636
    },
    "sat_max": {
      "Date": "2024-05-11T20:50:55",
      "MJD": 60441.87207,
      "ra": 140.1447,
      "dec": 71.225,
      "az": 334.2671,
      "alt": 53.201,
      "Distance": 513.29,
      "Separation": 36.799,
      "PA": 25,
      "Speed": 4.072
    }
  },
  "sunlit_duration_minutes": 10.23,
  "run_command": "sat_visibility stations.txt -T -l44.52324,11.33914,23 -d2024-05-11T20:00:00,2024-05-11T22:00:00 -s 25544",
  "status": 0,
  "errmsg": "",
  "n_sats_found": 1
}


  Where:
    sat_start/end/mid refer to times of Start/End/Mid of the visivility range,
    "max" for the highest satellite altitude over the horizon

    Distance: distance to satellite in km
    Separation: angular separation in degrees from the zenith
    PA: position angle of motion
    Speed: apparent angular rate of motion in arcminutes/second (== degrees/minute)


  LN @ INAF-OAS, Jan 2020.  Last change: 15/05/2024
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
#include "sat_visibility_def.h"


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
  "  -a Alt_min,Alt_max Geodetic altitude range filter (km)\n"
  "  -d CalendarDates   Calendar date (UTC) range of interest (in the form yyyy-mm-ddThh:mm:ss[.sss],yyyy-...)\n"
  "  -D DirTLEs         Directory with the repository of TLE files (def. ./; ignore -T option)\n"
  "  -i SatIntnlName    Satellite selection via its international designator\n"
  "  -j MJDs            Modified Julian Date range of interest (MJDstart,MJDend - ignored if Calendar Date given)\n"
  "  -l Lat,Lon,Alt     Geodetic observing site (comma separated data with no spaces)\n"
  "  -r radius          Region search radius - centred at the site zenith (degrees; def. 90)\n"
  "  -s SatNorad_n      Satellite selection via its NORAD number\n"
  "  -t deltaT          Visibility check delta time (seconds; def. 3)\n"

  "\nSwitches:\n"
  "  -h                 Print this help\n"
  "  -N                 Skip visibility check until Sun is set - if it is up at start date\n"
  "  -T                 Use default repository directory for TLE files (def. ./; see sat_visibility_def.h)\n\n"

  "Example usage:\n"
  "  %s stations.txt -T -l44.52324,11.33914,23 -d2024-05-11T20:00:00,2024-05-11T22:00:00 -s 25544\n"
  "  %s stations.txt -T -t 10 -l-29.25627,-70.73805,2400 -j60436.54,60436.55 -s 25544\n", progname, progname, progname);
}


void open_json() {
  printf("{"
    "\"swinfo\": {"
    "\"name\": \"%s\", \"author\": \"%s\", \"date\": \"%s\", \"version\": \"%s\"}, ",
    progname, progauthor, progdate, progversion);
}

void add_params_json() {
  printf("\"input_params\": {"
    "\"tle_file\": \"%s\", \"location\": {\"lat\":%8.4lf, \"lon\":%9.4lf, \"alt\":%8.1lf}, \"region\": {\"ra\":%9.4lf, \"dec\":%8.4lf, \"radius\":%8.4lf, \"lmst\":%8.4lf, \"az\":%9.4lf, \"alt\":%8.4lf, \"parang\":%8.3lf}, "
    "\"mjd_range\": [%11.5lf, %11.5lf], \"date_utc_range\": [\"%s\", \"%s\"], \"gmst\":%8.4lf, \"delta_time_s\": %d, "
    "\"notes\": \"All coordinates and radius in degrees. GMST, LMST in hrs.\"}, ",
    p.tle_file_name, p.lat, p.lon, p.ht_in_meters, p.ra_deg, p.de_deg, p.search_radius, p.lmst, p.az, p.alt, p.parang, p.mjd[0], p.mjd[1], p.date[0], p.date[1], p.gmst, p.delta_time);
}

void add_fieldsdesc_json() {
  printf(
    "\"sat_fields\": {\"desc\": [\"Date start/end sunlit\", \"MJD start/end sunlit\", \"RA J2000\", \"Dec J2000\", \"Azimuth\", \"Altitude\", \"Distance to sat.\", \"Angular separation\", \"Position angle\", \"Apparent angular rate of motion\"], "
    "\"type\": [\"datetime\", \"double\", \"double\", \"double\", \"double\", \"double\", \"double\", \"float\", \"float\", \"float\"], "
    "\"unit\": [\"\", \"MJD\", \"deg\", \"deg\", \"deg\", \"deg\", \"km\", \"deg\", \"deg\", \"arcmin/s\"]},");
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

void add_satlatlon_json(char *label) {
  printf("\"%s\": {\"lat\":%7.3lf, \"lon\":%7.3lf, \"alt\":%9.2lf, \"theta\":%8.3lf}, ",
    label,  geo.lat, geo.lon, geo.alt, geo.theta);
}

void add_satdata_json(char *label) {
  printf("\"%s\": {"
    "\"Date\": \"%s\", \"MJD\": %11.5lf, \"ra\": %8.4lf, \"dec\": %8.4lf, \"az\":%9.4lf, \"alt\":%8.4lf, \"Distance\": %9.2lf, \"Separation\": %6.4lf, \"PA\": %3d, \"Speed\": %6.3lf}",
    label,
    sat.cur_date, sat.cur_mjd, sat.ra * RAD2DEG, sat.dec * RAD2DEG, sat.az, sat.alt,
    sat.d_km, sat.ang_sep, (int)(sat.parang * RAD2DEG), sat.speed);
    if ( strstr(label,"sat_max") == NULL )
	printf(", ");
}

void add_sundata_json(char *label) {
  printf("\"%s\": {"
    "\"ra\":%7.3lf, \"dec\":%7.3lf, \"az\":%8.3lf, \"alt\":%7.3lf, \"lon\":%8.3lf, \"parang\":%8.3lf, \"separation_deg\":%7.3lf}, ",
    label,
    sun.ra, sun.dec, sun.az, sun.alt, sun.lon, sun.parang, sun.ang_sep);
}

void close_stat_json(int argc, char **argv, bool input_list, int status, char *errmsg, int n_sats_found) {
  int iarg = 1;
  printf("\"run_command\": \"%s", progname);
  if ( input_list ) {
	printf (" @%s", argv[1]);
	iarg++;
  }

  for ( int i = iarg; i < argc; i++ ) {
	printf (" %s", argv[i]);
  }
  printf("\", \"status\": %d, \"errmsg\": \"%s\", \"n_sats_found\": %d"
	 "}\n", status, errmsg, n_sats_found);
}


/*
  Compute intesection points distance on Earth of the Satellite --> Sun connecting line.
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



/*
  Compute satellite geodetic Lon, Lat, Alt and theta (RA) position given its ECI position and GMST.
  Returned values are deg, deg, km, deg.
*/

void sat_geoLocation(double *pos) {
  double r, e2, phi, sinphi, c;

  geo.theta = atan2(pos[1], pos[0]) * RAD2DEG;  /* degrees */
  geo.lon = (geo.theta - p.gmst * 15.);
  if (geo.lon > 180)
	geo.lon -= 360.;
  else if (geo.lon < -180.)
	geo.lon += 360.;

  r = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  e2 = F * (2. - F);
  geo.lat = atan2(pos[2], r);

  do {
	phi = geo.lat;
	sinphi = sin(phi);
	c = 1. / sqrt(1. - e2 * sinphi * sinphi);
	geo.lat = atan2(pos[2] + ERAD * c * e2 * sinphi, r);
  } while (fabs(geo.lat - phi) >= 1e-10);

  geo.alt = r/cos(geo.lat) - ERAD * c;  /* kilometers */

  geo.lat *= RAD2DEG;  /* degrees */
   if (geo.lat > 90.)
	geo.lat -= 360.;  /* could keep in -90, 90 range, but then need to manage geo->lon accordingly */
}


/* Compute the main Sun data and angular separation from given target at a given JD time */

void get_sun_data(double t, double target_ra_rad, double target_dec_rad, bool get_distance) {
  lpsun_radec(t, &sun.ra, &sun.dec);
/* Lon in [0, 360[ deg +W  (GEO_LON = -GHA = -(GMST - RA)) */
  sun.lon = fmod(sun.ra * RAD2DEG - p.gmst * 15. + 360., 360.);
  if ( sun.lon > 180. )  /* Lon in [-180, +180] deg +E */
	sun.lon = sun.lon - 360.;

  sun.ha = p.lmst - sun.ra * RAD2HRS;
  dechalat2alt(sun.dec, sun.ha, p.lat, &sun.alt, &sun.az, &sun.parang);
  sun.ang_sep = skysep_h(target_ra_rad, target_dec_rad, sun.ra, sun.dec);

  if ( get_distance ) {
/* Sun cartesian coordinates (in km) */ 
	sphrad2v(sun.ra, sun.dec, sun.d_km);

	for (int i = 0; i < 3; i++)
		sun.d_km[i] *= AU_KM;
  }

  sun.ra *= RAD2DEG;
  sun.dec *= RAD2DEG;
}


/* Compute posn_ang_of_motion and speed given initial coords and a time > t_initial */

void get_sat_PA_and_speed(double t_since) {
	double pos[3], ra1, dec1, d_ra, d_dec, dist_to_satellite;
	if ( sat.is_deep ) {
	  SDP4_init(sat_params, &tle);
	  SDP4(t_since, &tle, sat_params, pos, NULL);
	} else {
	  SGP4_init(sat_params, &tle);
	  SGP4(t_since, &tle, sat_params, pos, NULL);
	}
	get_satellite_ra_dec_delta(geo.observer_loc, pos, &ra1, &dec1, &dist_to_satellite);  /* Returned RA, Dec in radians */
	//-- epoch_of_date_to_j2000(t, &ra1, &dec1);  /* Do not apply precession */

	d_ra = (ra1 - sat.ra + TWOPI);
	if ( d_ra > PI )
		d_ra -= TWOPI;
	d_ra *= cos(sat.dec);
	d_dec = dec1 - sat.dec;
	sat.parang = atan2(d_ra, d_dec);
	if ( sat.parang < 0. )
		sat.parang += TWOPI;
	sat.speed = skysep_h(sat.ra, sat.dec, ra1, dec1) * 60 / p.delta_time;  /* arcminutes/second */
}



/* -- Main -- */

int main(int argc, char **argv)
{
  char errmsg[128] = {""}, tle_path[128] = {"."}, line1[100], line2[100],
	sat_name[25], tle_path_file[200], opt, *endptr;
  double jd[2], rho_sin_phi, rho_cos_phi,
	 target_ra_rad, target_dec_rad;
  int i, par_pos, status = 0, n_sats_found = 0,
	 len_intl_desig = 0;
  bool in_geopos = false,  /* make sure Geo location is passed */
	input_list = false;  /* if the input file is in fact a TLEs list */


  if ( argc < 2 ) {
	open_json();
	sprintf(errmsg, "Example usage:\n %s stations.txt -l-29.25627,-70.73805,2400 -j58941.72,58941.73 -s 25544", progname);
 	close_stat_json(argc, argv, input_list, 1, errmsg, n_sats_found);
	exit(0);
  }


/* Set default values */
  p.search_radius = 90.;  /* Default search radius in degrees */
  p.delta_time = 3;       /* Default time step (s) for second epoch pos. computation */
  p.skip_daytime = false;      /* If daytime must be ignored */
  p.single_sat_i = false;      /* If enquire for just 1 satellite using its international designator */
  p.single_sat_n = false;      /* If enquire for just 1 satellite using its NORAD number */
  p.use_deftledir = false;     /* If default directory with TLE files should be used (def. ./) */
  p.altrng_requested = false;  /* If altitude range filter requested */
  p.intl_desig[0] = '\0';
  p.satname[0] = '\0';
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
 	  case 'N':
		p.skip_daytime = true;
		i--;
		break;
 	  case 'T':
		p.use_deftledir = true;
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
		break;
          case 'd':
		if ( strchr(argv[i] + par_pos, '-') == NULL ) {  // A minimal check
			Usage();
		        exit(0);
		}
		sscanf(argv[i] + par_pos, "%19s,%19s", p.date[0], p.date[1]);
		date2mjd_nf(p.date[0], &p.mjd[0]);
		date2mjd_nf(p.date[1], &p.mjd[1]);
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
		break;
          case 'j':
		if ( strchr(argv[i] + par_pos, '-') != NULL ) {  // A minimal check
			Usage();
		        exit(0);
		}
		sscanf(argv[i] + par_pos, "%lf,%lf", &p.mjd[0], &p.mjd[1]);
		strcpy(p.date[0], mjd2date(p.mjd[0]));
		strcpy(p.date[1], mjd2date(p.mjd[1]));
		break;
          case 'l':
		sscanf(argv[i] + par_pos, "%lf,%lf,%lf", &p.lat, &p.lon, &p.ht_in_meters);
		/* Longitude locally in the range [0, 360[ (not [-180, +180]) degrees */
		if ( p.lon < 0. )
		  p.lon += 360;
		if ( p.lat < -90. || p.lat > 90. )  {
		  open_json();
		  sprintf(errmsg, "Latitude must be in the range [-90, +90] degrees. Read '%lf'", p.lat);
  		  close_stat_json(argc, argv, input_list, -1, errmsg, n_sats_found);
		  exit(-1);
		}
		in_geopos = true;
		break;
          case 'r':
		p.search_radius = atof(argv[i] + par_pos);
		break;
          case 's':
		p.norad_n = strtol(argv[i] + par_pos, &endptr, 10);
		p.single_sat_n = true;
		p.altrng_requested = false;  // No altitude filter
		break;
          case 't':
		p.delta_time = atof(argv[i] + par_pos);
		break;
          default:
		open_json();
		sprintf(errmsg, "Unrecognized command-line option '%s'", argv[i]);
  		close_stat_json(argc, argv, input_list, -2, errmsg, n_sats_found);
		exit(-2);
		break;
	}  // end switch
  }  // end if - for

  jd[0] = JD_OFF + p.mjd[0];  // From MJD to JD
  jd[1] = JD_OFF + p.mjd[1];

/* Use default TLE repository dir. ? */
  if ( p.use_deftledir )
	strcpy(tle_path, DEF_TLEDIR);

  if ( !in_geopos ) {
	Usage();
	exit(0);
  }

/* Figure out where the observer really is in Cartesian coordinates of date */
  earth_lat_alt_to_parallax(p.lat * DEG2RAD, p.ht_in_meters, &rho_cos_phi, &rho_sin_phi);
  observer_cartesian_coords(jd[0], p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, geo.observer_loc);
 
  const short fname_len = 64;
  char row[fname_len], norad_name[7], tle_list_file[200];
  FILE *lisfile, *ifile = NULL;

  open_json();


/* TLE file passed as first or last parameter? */
  int iarg = 1;
  if ( argv[1][0] == '-' )
	iarg = argc - 1;

  if ( *argv[iarg] == '@' ) {
	input_list = true;
	sprintf(tle_list_file, "%s/%s", tle_path, ++argv[iarg]);
	if ( !(lisfile = fopen(tle_list_file, "rb")) ) {
		sprintf(errmsg, "Could not open input file %s", tle_list_file);
  		close_stat_json(argc, argv, input_list, 2, errmsg, n_sats_found);
		exit(2);
	}
  }

  p.tle_file_name = argv[iarg];

  if ( !p.single_sat_n && !p.single_sat_i ) {
	sprintf(errmsg, "Need to input a NORAD ID or an international designator.");
  	close_stat_json(argc, argv, input_list, 2, errmsg, n_sats_found);
	exit(2);
  }

  double hazenith = 0.;
  char *tle_file_name;
  bool single_sat_found = false, read_tle_list = true;

  char intl_desig[12] = {"none"};

/* Loop on list of tle files or simple open of single input file */
  while ( read_tle_list ) {

    char *lastchar;
    if ( !input_list ) {
	if ( strchr(p.tle_file_name, '/') == NULL )
	  sprintf(tle_path_file, "%s/%s", tle_path, p.tle_file_name);
	else
	  strcpy(tle_path_file, p.tle_file_name);
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
  	close_stat_json(argc, argv, input_list, 2, errmsg, n_sats_found);
	exit(2);
    }

    if ( !fgets( line1, sizeof(line1), ifile) ) {
	sprintf(errmsg, "Could not read first TLE line from file %s", tle_path_file);
	close_stat_json(argc, argv, input_list, -2, errmsg, n_sats_found);
	exit(-2);
    }

    strncpy(sat_name, line1, 23);
    trimend(sat_name);

    while ( fgets(line2, sizeof(line2), ifile) ) {

      if ( !parse_elements(line1, line2, &tle) ) { /* TLE found */
/* First check if we have processed this sat already */
	strncpy(norad_name, line1+2, 5);  /* this is tle.norad_number */
	norad_name[5] = '\0';

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
	}
 
	if ( single_sat_found ) {
		break;
	}

      }  // end TLE found

      strcpy(line1, line2);
      if ( line2[0] != '1' && line2[0] != '2' ) {
	strncpy(sat_name, line2, 23);
	trimend(sat_name);
      }
    }  // end while fgets(line2 ...

    if ( single_sat_found )
	break;

  }  // end while read_tle_list

  fclose(ifile);

  if ( !single_sat_found ) {
	sprintf(errmsg, "Satellite not found in %s", tle_path_file);
	close_stat_json(argc, argv, input_list, 1, errmsg, n_sats_found);
	exit(1);
  }


   sat.is_deep = select_ephemeris(&tle);

  double t_since, t_min_ang_dist = 0, min_ang_sep = 181, parang_sat,
	mjd_sunlit_start = 0, mjd_sunlit_mid = 0, mjd_sunlit_end = 0,
	pos[3],  /* Satellite position vector */
	t = jd[0],
	eray[2];  /* Sun-Earth ray intersection radius */

  bool is_t_start = true, sunlit_start_found = false, sunlit_end_found = false;


  p.de_deg = p.lat;
  target_dec_rad = p.lat * DEG2RAD;

  add_geofieldsdesc_json();
/* Alt, Az, Parang from Dec, HA, Lat */
  dechalat2alt(target_dec_rad, hazenith, p.lat, &p.alt, &p.az, &p.parang);
  add_params_json();


/* Recursive computation of Sun and satellite position at step of "delta_time" seconds */
/* TODO: compute Sun rise-set time to check for meaningful requested time range */


/* Skip daytime */
  if ( p.skip_daytime ) {
// Loop at 1 min step until Sun is set
    while (1) {
	p.lmst = lmst_hr(t, p.lon, &p.gmst);
	lpsun_radec(t, &sun.ra, &sun.dec);
	sun.ha = p.lmst - sun.ra * RAD2HRS;
	dechalat2alt(sun.dec, sun.ha, p.lat, &sun.alt, &sun.az, &sun.parang);
	if ( sun.alt < 0. || t > jd[1] )
		break; 
	t += 0.0007;
    }
  }

  while (t <= jd[1]) {  // Final t is approx

/* Local Mean Sidereal Time & GMST */
    p.lmst = lmst_hr(t, p.lon, &p.gmst);

/* Always use zenith coords */
    p.ra_deg = p.lmst * 15.;  
    target_ra_rad = p.lmst / RAD2HRS;

    t_since = (t - tle.epoch) * 1440.;  /* From days to minutes */

    if ( sat.is_deep ) {
	SDP4_init(sat_params, &tle);
	SDP4(t_since, &tle, sat_params, pos, NULL);
    } else {
	SGP4_init(sat_params, &tle);
	SGP4(t_since, &tle, sat_params, pos, NULL);
    }

    if ( p.altrng_requested ) {
	sat_geoLocation(pos);

/* If altitude filter requested */
	if ( geo.alt < p.alt_min || geo.alt > p.alt_max ) {
	  sprintf(errmsg, "Satellite altitude out of requested range");
	  close_stat_json(argc, argv, input_list, 2, errmsg, n_sats_found);
	  exit(1);
	}
    }

    observer_cartesian_coords(t, p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, geo.observer_loc);

    get_satellite_ra_dec_delta(geo.observer_loc, pos, &sat.ra, &sat.dec, &sat.d_km);

    sat.ang_sep = skysep_h(target_ra_rad, target_dec_rad, sat.ra, sat.dec);

/* Compute initial speed to be used to accelerate visibility search (TODO) */
    if ( is_t_start ) {
	is_t_start = false;
	get_sat_PA_and_speed(t_since);
    }

/*
printf("\nMJD, p.lmst, t_since, tle.epoch: %lf, %lf, %lf, %lf\n",t-JD_OFF, p.lmst, t_since, tle.epoch);
printf("target_ra_rad, target_dec_rad, sat.ra, sat.dec, sat.ang_sep: %lf, %lf, %lf, %lf,  %lf\n", target_ra_rad*RAD2DEG, target_dec_rad*RAD2DEG, sat.ra*RAD2DEG, sat.dec*RAD2DEG, sat.ang_sep);
double ra_tmp=ra, de_tmp=dec;
epoch_of_date_to_j2000(t, &ra_tmp, &de_tmp);
printf("raj2000, decj2000: %lf, %lf\n", ra_tmp*RAD2DEG, de_tmp*RAD2DEG);
*/

    if ( sat.ang_sep > p.search_radius ) {
	t += p.delta_time / SEC_IN_DAY;  // add deltaT
	if ( !sunlit_start_found )
	  continue;
	break;
    }

    get_sun_data(t, target_ra_rad, target_dec_rad, true);

/* If satellite is sunlit (TODO to account for Sun angular size) */
    sat.is_sunlit = ! intersect_satsun_sphere(pos, sun.d_km, eray);
//printf("\n%d - %s to Sun ray intersect Earth at %lf  %lf\n", (int)sat.is_sunlit, sat_name, eray[0], eray[1]);

    if ( sat.ang_sep < min_ang_sep && sat.is_sunlit ) {
	min_ang_sep = sat.ang_sep;
	t_min_ang_dist = t;
    }

    t += p.delta_time / SEC_IN_DAY;  // add deltaT

    if ( !sat.is_sunlit ) {
	if ( !sunlit_start_found )
	  continue;

/* First time the satellite no longer sunlit */
	if ( !sunlit_end_found ) {
	  sunlit_end_found = true;
	  mjd_sunlit_end = t - JD_OFF;
	  sat_geoLocation(pos);
	  add_satlatlon_json("geoloc_end");
	  add_sundata_json("sun_end");
	  sat.cur_mjd = t - JD_OFF - p.delta_time / SEC_IN_DAY;
	  strcpy(sat.cur_date, mjd2date(sat.cur_mjd));
	  sat.ha = p.lmst - sat.ra * RAD2HRS;
	  dechalat2alt(sat.dec, sat.ha, p.lat, &sat.alt, &sat.az, &parang_sat);
	  epoch_of_date_to_j2000(t, &sat.ra, &sat.dec);  /* Approx precession */
	  add_satdata_json("sat_end");
	  break;
      }
    }

/* First time the satellite is sunlit */
    if ( !sunlit_start_found ) {
	sunlit_start_found = true;
	mjd_sunlit_start = t - JD_OFF;

	if ( !p.single_sat_i ) {
		if ( tle.intl_desig[0] != ' ' )
		  snprintf( intl_desig, 10, "%s%.2s-%s",
			(atoi(tle.intl_desig) > 57000 ? "19" : "20"),
			tle.intl_desig, tle.intl_desig + 2);
	}

	t_since += p.delta_time * 1440. / SEC_IN_DAY;

	get_sat_PA_and_speed(t_since);

	add_fieldsdesc_json();
	printf(" \"satellite\": ");
	printf("{\"name\": \"%s\", \"intl_desig\": \"%s\", \"norad_n\": %d, ",
		sat_name, intl_desig, tle.norad_number);

	sat_geoLocation(pos);
	add_satlatlon_json("geoloc_start");
	add_sundata_json("sun_start");
	sat.cur_mjd = t - JD_OFF - p.delta_time / SEC_IN_DAY;
	strcpy(sat.cur_date, mjd2date(sat.cur_mjd));
	sat.ha = p.lmst - sat.ra * RAD2HRS;
	dechalat2alt(sat.dec, sat.ha, p.lat, &sat.alt, &sat.az, &parang_sat);
	epoch_of_date_to_j2000(t, &sat.ra, &sat.dec);  /* Approx precession */
	add_satdata_json("sat_start");
	  n_sats_found++;

    } // end if !sunlit_start_found

  }  // end time loop


  t -= p.delta_time / SEC_IN_DAY;

/* Still sunlit at the end of the requested time range */
  if ( sunlit_start_found && !sunlit_end_found ) {
	mjd_sunlit_end = t - JD_OFF;
	sat_geoLocation(pos);
	add_satlatlon_json("geoloc_end");
	add_sundata_json("sun_end");
	sat.cur_mjd = mjd_sunlit_end;
	strcpy(sat.cur_date, mjd2date(sat.cur_mjd));
	sat.ha = p.lmst - sat.ra * RAD2HRS;
	dechalat2alt(sat.dec, sat.ha, p.lat, &sat.alt, &sat.az, &parang_sat);
	epoch_of_date_to_j2000(t, &sat.ra, &sat.dec);  /* Approx precession */
	add_satdata_json("sat_end");
  }


/* Compute data at Tmid */
  if ( sunlit_start_found ) {
	mjd_sunlit_mid = mjd_sunlit_start + (mjd_sunlit_end - mjd_sunlit_start)/2.;
	t = mjd_sunlit_mid + JD_OFF;
	p.lmst = lmst_hr(t, p.lon, &p.gmst);
	p.ra_deg = p.lmst * 15.;
	target_ra_rad = p.lmst / RAD2HRS;

	t_since = (t - tle.epoch) * 1440.;
	
	if ( sat.is_deep ) {
	  SDP4_init(sat_params, &tle);
	  SDP4(t_since, &tle, sat_params, pos, NULL);
	} else {
	  SGP4_init(sat_params, &tle);
	  SGP4(t_since, &tle, sat_params, pos, NULL);
	}

	sat_geoLocation(pos);

	observer_cartesian_coords(t, p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, geo.observer_loc);

	get_satellite_ra_dec_delta(geo.observer_loc, pos, &sat.ra, &sat.dec, &sat.d_km);

	sat.ang_sep = skysep_h(target_ra_rad, target_dec_rad, sat.ra, sat.dec);

	get_sun_data(t, target_ra_rad, target_dec_rad, false);
	add_satlatlon_json("geoloc_mid");
	add_sundata_json("sun_mid");
	strcpy(sat.cur_date, mjd2date(mjd_sunlit_mid));
	sat.ha = p.lmst - sat.ra * RAD2HRS;
	dechalat2alt(sat.dec, sat.ha, p.lat, &sat.alt, &sat.az, &parang_sat);
	epoch_of_date_to_j2000(t, &sat.ra, &sat.dec);  /* Approx precession */
	add_satdata_json("sat_mid");
//printf("\n\nMid: %lf, %lf, %lf, %lf\n", sat.ha, sat.alt, sat.az, parang_sat);


/* Compute data at minimun Zenith distance (maximum altitude) */
	t = t_min_ang_dist;

	p.lmst = lmst_hr(t, p.lon, &p.gmst);
	p.ra_deg = p.lmst * 15.;
	target_ra_rad = p.lmst / RAD2HRS;

	t_since = (t - tle.epoch) * 1440.;
	
	if ( sat.is_deep ) {
	  SDP4_init(sat_params, &tle);
	  SDP4(t_since, &tle, sat_params, pos, NULL);
	} else {
	  SGP4_init(sat_params, &tle);
	  SGP4(t_since, &tle, sat_params, pos, NULL);
	}

	sat_geoLocation(pos);

	observer_cartesian_coords(t, p.lon * DEG2RAD, rho_cos_phi, rho_sin_phi, geo.observer_loc);

	get_satellite_ra_dec_delta(geo.observer_loc, pos, &sat.ra, &sat.dec, &sat.d_km);

	sat.ang_sep = skysep_h(target_ra_rad, target_dec_rad, sat.ra, sat.dec);

	get_sun_data(t, target_ra_rad, target_dec_rad, false);
	add_satlatlon_json("geoloc_max");
	add_sundata_json("sun_max");
	strcpy(sat.cur_date, mjd2date(t - JD_OFF));
	sat.ha = p.lmst - sat.ra * RAD2HRS;
	dechalat2alt(sat.dec, sat.ha, p.lat, &sat.alt, &sat.az, &parang_sat);
	epoch_of_date_to_j2000(t, &sat.ra, &sat.dec);  /* Approx precession */
	add_satdata_json("sat_max");

  }


  if ( n_sats_found > 0 ) {
	printf("},"
	" \"sunlit_duration_minutes\": %.2f, ", (float)(mjd_sunlit_end - mjd_sunlit_start)*1440.);
  } else {
	p.lmst = lmst_hr(jd[0], p.lon, &p.gmst);
	target_ra_rad = p.lmst / RAD2HRS;
	get_sun_data(jd[0], target_ra_rad, target_dec_rad, false);
	add_sundata_json("sun");
  }


  close_stat_json(argc, argv, input_list, status, errmsg, n_sats_found);

  return(0);
}
