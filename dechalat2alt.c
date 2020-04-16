/*
   Returns altitude (deg) for Dec, HA, Lat (rad, hr, deg).
   Also computes and returns azimuth through pointer argument,
   and as an extra added bonus returns parallactic angle (deg)
   through another pointer argument.
 
*/

#include <math.h>

#include "const_def.h"

double dechalat2alt(double dec, double ha, double lat, double *az, double *parang)
{
        double x,y,z;
        double sinp, cosp;  /* sin and cos of parallactic angle */
        double cosdec, sindec, cosha, sinha, coslat, sinlat;

        ha = ha * HRS2RAD;
        lat = lat * DEG2RAD;  /* thank heavens for pass-by-value */
        cosdec = cos(dec); sindec = sin(dec);
        cosha = cos(ha); sinha = sin(ha);
        coslat = cos(lat); sinlat = sin(lat);
        x = RAD2DEG * asin(cosdec*cosha*coslat + sindec*sinlat);
        y =  sindec*coslat - cosdec*cosha*sinlat; /* due N comp. */
        z =  -1. * cosdec*sinha; /* due east comp. */
        *az = atan2(z,y);

        /* as it turns out, having knowledge of the altitude and
           azimuth makes the spherical trig of the parallactic angle
           less ambiguous ... so do it here!  Method uses the
           "astronomical triangle" connecting celestial pole, object,
           and zenith ... now know all the other sides and angles,
           so we can crush it ... */

        if(cosdec != 0.) { /* protect divide by zero ... */
           sinp = -1. * sin(*az) * coslat / cosdec;
                /* spherical law of sines .. note cosdec = sin of codec,
 *                         coslat = sin of colat .... */
           cosp = -1. * cos(*az) * cosha - sin(*az) * sinha * sinlat;
                /* spherical law of cosines ... also transformed to local
 *                       available variables. */
           *parang = atan2(sinp,cosp) * RAD2DEG;
                /* let the library function find the quadrant ... */
        }
        else { /* you're on the pole */
           if(lat >= 0.) *parang = 180.;
           else *parang = 0.;
        }

        *az *= RAD2DEG;  /* done with taking trig functions of it ... */
        while(*az < 0.) *az += 360.;  /* force 0 -> 360 */
        while(*az >= 360.) *az -= 360.;

        return(x);
}
