/*
   Definitions for sat_skymap
  
   LN @ INAF-OAS, Feb 2020                      Last change: 24/04/2020
*/

#ifndef SAT_SKYMAP_DEF_H
#define SAT_SKYMAP_DEF_H

const char *progname = "sat_skymap", *progauthor = "L. Nicastro @ INAF-OAS", *progdate = "2020-04-24", *progversion = "0.2b";

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
