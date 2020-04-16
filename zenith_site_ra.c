/*
   Equatorial RA coordinate (J2000) of the zenith at a given site and date.

   LN adapted
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "const_def.h"

void earthtilt (double tjd,double *mobl, double *tobl, double *eq,
  double *dpsi, double *deps);

void ut2gmst(double ut1, double *gmst)
{
  double jd0;    /* Julian day at midnight Universal Time */
  double secs;   /* Time of day, UT seconds since UT midnight */
  double T0, msday;

/* Julian day at given UT */
  jd0 = floor(ut1);
  secs = ut1 - jd0;
  if ( secs < 0.5 ) {
	jd0 -= 0.5;
	secs += 0.5;
   } else {
	jd0 += 0.5;
	secs -= 0.5;
  }
  secs *= SEC_IN_DAY;

  T0 = (jd0 - JD2000)/36525.0;

  /* Greenwich Mean Sidereal Time at 0h UT of date */
  *gmst = (((-2.0e-6*T0 - 3.e-7)*T0 + 9.27701e-2)*T0 + 8640184.7942063F)*T0
          + 24110.54841F;
#ifdef DEBUG
printf("GMST0: %lf\n", *gmst);
#endif
  msday = (((-(4. * 2.0e-6)*T0 - (3. * 3.e-7))*T0 + (2. * 9.27701e-2))*T0
          + 8640184.7942063F)/(SEC_IN_DAY*36525.) + 1.0;

/* GMST at this UT (in hours) */
  *gmst += msday*secs;
  *gmst = *gmst - SEC_IN_DAY * floor(*gmst/SEC_IN_DAY);

  *gmst /= 3600.;

  if ( *gmst < 0. )
	*gmst += 24.;
  else if ( *gmst >= 24. )
	*gmst -= 24.;
}


void ut2gast(double ut1, double *gast)
{
  double gmst, a, b, c, d, eqeq;

  ut2gmst(ut1, &gmst);

// Quantities related to the orientation of the Earth's rotation axis at 'jd'
  earthtilt (ut1, &a, &b, &eqeq, &c, &d);

// GAST at this UT (in hours)
  *gast = gmst + eqeq/3600.;

#ifdef DEBUG
printf("gmst, eqeq, gast: %lf, %lf, %lf\n", gmst, eqeq/3600., *gast);
#endif
}


void ut2last(double ut1, double lon, double *last)
{
  double gast;

  ut2gast(ut1, &gast);

//printf("GAST=%f\n",gast);
// LAST at this UT (in hours)
  *last = gast + lon/15.;

  if ( *last < 0. )
	*last += 24.;
  else if ( *last >= 24. )
	*last -= 24.;
}


// Input UT1 in MJD, Longitude in hours.
// Return RA in degrees.
int zenith_site_ra(double ut1, double lon, double *radeg)
{
  double gmst, lmst;

  ut1 += JD_OFF;  // From MJD to JD

  ut2gmst(ut1, &gmst);
  lmst = gmst + lon/15.;
#ifdef DEBUG
printf("gmst, lmst: %lf, %lf\n", gmst, lmst);
#endif

  if ( lmst < 0. )
	lmst += 24.;
  else if ( lmst >= 24. )
	lmst -= 24.;

  *radeg = lmst * 15.;  // - har; /* degrees */

  return 0;
}

// Version using LAST instead of LST.
int zenith_site_ra_last(double ut1, double lon, double *radeg)
{
  double last;

  ut1 += JD_OFF;  // From MJD to JD

  ut2last(ut1, lon, &last); /* From Universal Time to Local Apparent Sideral Time */

  *radeg = last * 15.;  // - har;

  return 0;
}


//
// Unused
void slaDh2e ( double az, double el, double phi, double *ha, double *dec )
{
  double sa, ca, se, ce, sp, cp, x, y, z, r;

/* Useful trig functions */
  sa = sin ( az );
  ca = cos ( az );
  se = sin ( el );
  ce = cos ( el );
  sp = sin ( phi );
  cp = cos ( phi );

/* HA,Dec as x,y,z */
  x = - ca * ce * sp + se * cp;
  y = - sa * ce;
  z = ca * ce * cp + se * sp;

/* To spherical */
  r = sqrt ( x * x + y * y );
  *ha = ( r == 0.0 ) ? 0.0 : atan2 ( y, x ) ;
  *dec = atan2 ( z, r );
}
