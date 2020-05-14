/*
   Low precision formulae for the sun, from Almanac p. C24 (1990).
   RA and Dec are returned in radians.
*/

#include <math.h>
  
#include "const_def.h"

/* Returns radian angle (0 - 2pi) for coords x, y (get that quadrant right !) */
double atan_circ(double x, double y)
{
  double theta;

  if ( x == 0. ) {
	if ( y > 0. )
       	  theta = M_PI_2;
	else if ( y < 0. )
       	  theta = 3.* M_PI_2;
	else
       	  theta = 0.;  /* x and y zero */
  }
  else theta = atan(y/x);

  if ( x < 0. )
	theta += PI;
  if ( theta < 0. )
	theta += TWOPI;

  return(theta);
}


void lpsun_radec(double jd, double *ra, double *dec)
{
  double n, g, lambda,epsilon,x,y,z;

  n = jd - JD2000;
  g = (357.528 + 0.9856003 * n) * DEG2RAD;
  lambda = (280.460 + 0.9856474 * n + 1.915 * sin(g) + 0.020 * sin(2. * g)) * DEG2RAD;
  epsilon = (23.439 - 0.0000004 * n) * DEG2RAD;

  x = cos(lambda);
  y = cos(epsilon) * sin(lambda);
  z = sin(epsilon) * sin(lambda);

  *ra = atan_circ(x,y);
  *dec = asin(z);
}
