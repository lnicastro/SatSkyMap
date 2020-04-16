/*
   Returns the local MEAN sidereal time (dec hrs) at julian date jd
   --- at west longitude
   +++ at Est longitude
   (decimal degrees).
   Follows definitions in 1992 Astronomical Almanac, pp. B7 and L2.
   Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
   p.359, 1982.  On workstations, accuracy (numerical only!)
   is about a millisecond in the 1990s.

*/

#include <math.h>

#include "const_def.h"

double lst_hr(const double jd, double longit)
{
	double t, ut, jdmid, jdint, jdfrac, sid_g;
        long jdin, sid_int;

        jdin = jd;         /* fossil code from earlier package which
                        split jd into integer and fractional parts ... */
        jdint = jdin;
        jdfrac = jd - jdint;
        if(jdfrac < 0.5) {
                jdmid = jdint - 0.5;
                ut = jdfrac + 0.5;
        }
        else {
                jdmid = jdint + 0.5;
                ut = jdfrac - 0.5;
        }
        t = (jdmid - JD2000)/36525;
        sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/SEC_IN_DAY;
        sid_int = sid_g;
        sid_g = sid_g - (double) sid_int;
        //sid_g = sid_g + 1.0027379093 * ut - longit/24.;
        sid_g = sid_g + 1.0027379093 * ut + longit/360.;
        sid_int = sid_g;
        sid_g = (sid_g - (double) sid_int) * 24.;
        if(sid_g < 0.)
		sid_g += 24.;
        return(sid_g);
}
