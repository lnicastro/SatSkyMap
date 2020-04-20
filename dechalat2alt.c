/*
   Returns altitude (deg) for Dec, HA, Lat (rad, hr, deg).
   Also computes and returns azimuth through pointer argument,
   and as an extra added bonus returns parallactic angle (deg)
   through another pointer argument.
 
*/

#include <math.h>

#include "const_def.h"

void dechalat2alt(double dec, double ha, double lat, double *alt, double *az, double *parang)
{
  double cosha, sinha, sinlat;
  double sinp, cosp;  /* sin and cos of parallactic angle */
  double sinaz, cosaz;

  lat = lat * DEG2RAD;
  double cosdec = cos(dec), sindec = sin(dec), coslat = cos(lat);


/* Manage peculiar cases and calculate Azimuth */
  if ( ha == 0. ) {
	*az = 0.;
	cosha = 1.;
	cosaz = 1.;
	sinha = 0.;
	sinaz = 0.;
	sinp  = 0.;
	cosp = -1;
  } else {
	ha = ha * 15 * DEG2RAD;
	cosha = cos(ha);
	sinha = sin(ha);
	if (coslat < 1e-6) {  /* very close to a pole: assume +/- 90 deg */
	  coslat = 0.;
	  if (lat < 0.)
	    sinlat = -1.;
	  else
	    sinlat = 1.;
          *az = ha;
	  sinaz = sinha;
	  cosaz = cosha;
	  sinp  = 0.;
	  cosp  = -1. * cosha*cosha - sinha*sinha;
	} else {
	  sinlat = sin(lat);

          *az = atan2(-1. * cosdec*sinha, sindec*coslat - cosdec*cosha*sinlat);  /* due east, due north comp. */

	  sinaz = sin(*az);
	  cosaz = cos(*az);

	  *az *= RAD2DEG;  /* done with taking trig functions of it ... */
          if (*az < 0.) *az += 360.;  /* force 0 -> 360 */
          if (*az >= 360.) *az -= 360.;
          if (cosdec != 0.) { /* protect divide by zero ... */
	    sinp = -1. * sinaz * coslat / cosdec;
	    cosp = -1. * cosaz * cosha - sinaz * sinha * sinlat;
	  }
	}
  }

/* Parallactic angle */
  if (cosdec != 0.) { /* protect divide by zero ... */
	*parang = atan2(sinp,cosp) * RAD2DEG;  /* let the library function find the quadrant */
  } else {  /* pole */
	if (lat >= 0.)
	   *parang = 180.;
	else
	   *parang = 0.;
   }

/* Altitude */
  if ( fabs(dec - lat) < 1e-6 )
	*alt = 90.;
  else
	*alt = RAD2DEG * asin(cosdec*cosha*coslat + sindec*sinlat);
}
