/*
   Definitions for sat_skymap
  
   LN @ INAF-OAS, Feb 2020                      Last change: 11/01/2023
*/

#ifndef SAT_SKYMAP_DEF_H
#define SAT_SKYMAP_DEF_H

/* Earth radius in km - IERS Conventions (2003) */
const double ERAD = 6.3781366e3;

/* Earth ellipsoid flattening factor - IERS Conventions (2003). Note: in "predict" it is 3.35281066474748E-3 */
const double F = 3.352819697896193e-3;

/* AU in km - IAU 2012 Resolution B2 */
const double AU_KM = 149597870.700; 

const char *progname = "sat_skymap",
	*progauthor = "L. Nicastro @ INAF-OAS",
       	*progdate = "2023-01-11",
       	*progversion = "0.3e";

  typedef struct myParams {
	  char *tle_file_name,
		date[24],
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
		mjd,
		alt_min,
		alt_max;
	    int delta_time,
		max_sats,
		norad_n;
	   bool haversine,
		info_only,
		sunlit_only,
		single_sat_i,
		single_sat_n,
		satname_filter,
		geoloc_requested,
		geoloc_reference,
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
		sep;
	double d_km[3];
  } Sun;

  typedef struct myGeoloc {
	 double lon,
		lat,
		alt,
		theta;
  } Geoloc;

/* Default TLE files repository */
#define DEF_TLEDIR "/usr/local/TLErepo"

#ifdef __cplusplus
extern "C" {
#endif

int date2mjd(char *date, double *mjd);
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
