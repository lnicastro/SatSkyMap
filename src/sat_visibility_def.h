/*
   Definitions for sat_visibility
  
   LN @ INAF-OAS, Feb 2020.  Last change: 13/05/2024
*/

#ifndef SAT_VISIBILITY_DEF_H
#define SAT_VISIBILITY_DEF_H

/* Earth radius in km - IERS Conventions (2003) */
const double ERAD = 6.3781366e3;

/* Earth ellipsoid flattening factor - IERS Conventions (2003). Note: in "predict" it is 3.35281066474748E-3 */
const double F = 3.352819697896193e-3;

/* AU in km - IAU 2012 Resolution B2 */
const double AU_KM = 149597870.700; 

const char *progname = "sat_visibility",
	*progauthor = "L. Nicastro @ INAF-OAS",
       	*progdate = "2024-05-13",
       	*progversion = "0.1a";

  typedef struct myParams {
	  char *tle_file_name,
		date[2][24],
		intl_desig[12],
		satname[24];
	 double lat,
		lon,
		ht_in_meters,
		ra_deg,
		de_deg,
		search_radius,
		az,
		alt,
		parang,
		gmst,
		lmst,
		mjd[2],
		alt_min,
		alt_max;
	    int delta_time,
		norad_n;
	   bool skip_daytime,
		single_sat_i,
		single_sat_n,
		altrng_requested,
		use_deftledir;
  } Params;

  typedef struct mySun {
	 double ra,
		dec,
		ha,
		az,
		alt,
		parang,
		lon,
		ang_sep,
		d_km[3];
  } Sun;

  typedef struct mySat {
	   char cur_date[24];
	double  cur_mjd,
		ra,
		dec,
		ha,
		az,
		alt,
		parang,
		speed,    /* Speed in arcminutes/second (== degrees/minute) */
		ang_sep,  /* Zenith distance in degrees */
		d_km;     /* Distance in km */
	   bool is_deep,  /* Deep space sat. */
		is_sunlit;
  } Sat;

  typedef struct myGeoloc {
	 double lon,
		lat,
		alt,
		theta,
		observer_loc[3];
  } Geoloc;


/* Initial / input parmaeters */
  Params p;
/* Satellite position */
  Sat sat;
/* Sun position */
  Sun sun;
/* Geodetic position */
  Geoloc geo;

  tle_t tle;  /* Structure for two-line elements set for satellite */
  double sat_params[N_SAT_PARAMS];

/* Default TLE files repository */
#define DEF_TLEDIR "/usr/local/TLErepo"

#ifdef __cplusplus
extern "C" {
#endif

int date2mjd_nf(char *date, double *mjd);
char *mjd2date(double mjd);

/* Haversine formula for spherical distance */
double skysep_h(const double theta1, const double phi1, const double theta2, const double phi2);
double lmst_hr(const double jd, double lon, double *gmst);
void dechalat2alt(double dec, double ha, double lat, double *alt, double *az, double *parang);
void lpsun_radec(double jd, double *ra, double *dec);

#ifdef __cplusplus
}
#endif

#endif
