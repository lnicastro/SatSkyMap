/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

/*
   LN adaped 
*/

#include <stdio.h>
#include <math.h>
#include "observe.h"

#include "const_def.h"

/*  Assorted functions useful in conjunction with the satellite code
   library for determining where the _observer_,  as well as the _target_,
   happens to be.  Combine the two positions,  and you can get the
   distance/RA/dec of the target as seen by the observer. */

#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)

/* function for Greenwich mean sidereal time, ripped from 'deep.cpp' */

static inline double ThetaG( double jd)
{
/* Reference:  The 1992 Astronomical Almanac, page B6. */
  const double omega_E = 1.00273790934;  /* Earth rotations per sidereal day (non-constant) */
  const double UT = fmod(jd + .5, 1.);
  double t_cen, GMST;

  t_cen = (jd - UT - JD2000) / 36525.;
  GMST = 24110.54841 + t_cen * (8640184.812866 + t_cen * (0.093104 - t_cen * 6.2E-6));
#ifdef DEBUG
printf("GMST0: %lf\n", GMST);
#endif

  GMST = fmod( GMST + SEC_IN_DAY * omega_E * UT, SEC_IN_DAY);
  if( GMST < 0.)
     GMST += SEC_IN_DAY;

#ifdef DEBUG
printf("ThetaG-gmst: %lf\n", 24. * GMST/SEC_IN_DAY);
#endif


  return(TWOPI * GMST / SEC_IN_DAY);
} /*Function thetag*/


void DLL_FUNC observer_cartesian_coords( const double jd, const double lon,
              const double rho_cos_phi, const double rho_sin_phi,
              double *vect)
{
   const double lmst = lon + ThetaG(jd);

   *vect++ = cos(lmst) * rho_cos_phi * EARTH_MAJOR_AXIS / 1000.;
   *vect++ = sin(lmst) * rho_cos_phi * EARTH_MAJOR_AXIS / 1000.;
   *vect++ = rho_sin_phi             * EARTH_MAJOR_AXIS / 1000.;
}


void DLL_FUNC earth_lat_alt_to_parallax( const double lat, const double ht_in_meters,
                    double *rho_cos_phi, double *rho_sin_phi)
{
   const double u = atan( sin( lat) * EARTH_AXIS_RATIO / cos( lat));

   *rho_sin_phi = EARTH_AXIS_RATIO * sin( u) +
                           (ht_in_meters / EARTH_MAJOR_AXIS) * sin( lat);
   *rho_cos_phi = cos( u) + (ht_in_meters / EARTH_MAJOR_AXIS) * cos( lat);
}


void DLL_FUNC get_satellite_ra_dec_delta( const double *observer_loc,
                                 const double *satellite_loc, double *ra,
                                 double *dec, double *delta)
{
   double vect[3], dist2 = 0.;
   int i;

   for( i = 0; i < 3; i++)
      {
      vect[i] = satellite_loc[i] - observer_loc[i];
      dist2 += vect[i] * vect[i];
      }
   *delta = sqrt( dist2);
   *ra = atan2( vect[1], vect[0]);
   if( *ra < 0.)
      *ra += PI + PI;
   *dec = asin( vect[2] / *delta);
}


/* Formulae from Meeus' _Astronomical Algorithms_ for approximate precession.
More than accurate enough for our purposes.  */

void DLL_FUNC epoch_of_date_to_j2000( const double jd, double *ra, double *dec)
{
  const double t_centuries = (jd - JD2000) / 36525.;
  const double m = (3.07496 + .00186 * t_centuries / 2.) * DEG2RAD / 240.;
  const double n = (1.33621 - .00057 * t_centuries / 2.) * DEG2RAD / 240.;
  const double ra_rate  = m + n * sin( *ra) * tan( *dec);
  const double dec_rate = n * cos( *ra);

  *ra -= t_centuries * ra_rate * 100.;
  *dec -= t_centuries * dec_rate * 100.;

/* LN add: RA in 0 - 2pi range */
  //*ra = fmod(*ra + PI*10., TWOPI);
  if ( *ra < 0. ) 
	*ra += TWOPI;
  else if ( *ra >= TWOPI )
	*ra -= TWOPI;
}
