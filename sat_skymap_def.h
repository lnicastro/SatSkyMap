/*
 * Definitions for sat_skymap
 *
 *  LN @ INAF-OAS, Feb 2020                      Last change: 19/04/2020
 *
 */

#ifndef SAT_SKYMAP_DEF_H
#define SAT_SKYMAP_DEF_H

const char *progname = "sat_skymap", *progauthor = "L. Nicastro @ INAF-OAS", *progdate = "2020-04-19", *progversion = "0.2a";

/* Default TLE files repository */
#define TLE_DATADIR "/usr/local/TLErepo/"

int date2mjd(char *date, double *mjd);
char *mjd2date(double mjd);

/* Haversine formula for spherical distance */
double skysep_h(const double theta1, const double phi1, const double theta2, const double phi2);
double lst_hr(const double jd, double lon, double *gmst);
void dechalat2alt(double dec, double ha, double lat, double *alt, double *az, double *parang);
void lpsun_radec(double jd, double *ra, double *dec);

#endif
