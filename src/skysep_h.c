#include <math.h>

#include "const_def.h"

/*
  Haversine formula for spherical distance  (R.W. Sinnot, 1984).

  Note: input angles in radians, returned separation in degrees.
*/
double skysep_h(const double theta1, const double phi1, const double theta2, const double phi2)
{
  double radif, sin2a, sin2d,
	 t1 = theta1, t2 = theta2;

// Manage Geodetic Lon which is in [-180, +180]. Do not check Phi
  if ( theta1 < 0. )
	t1 += TWOPI;
  if ( theta2 < 0. )
	t2 += TWOPI;

  radif  = fabs(theta2-theta1);
  if (radif > PI) radif = TWOPI - radif;
  sin2a = sin(radif/2.);
  sin2a *= sin2a;
  sin2d = sin((phi2-phi1)/2.);
  sin2d *= sin2d;
  return (2 * asin( sqrt(sin2d + cos(phi1)*cos(phi2)*sin2a) ) * RAD2DEG);
}

