/*
  From DATE in the form "2020-08-07T12:15:11" to MJD.

  Return 0 un success.


  LN@INAF-OAS, Aug 2006   Last change: 06/05/2024
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static const double JD_OFF = 2400000.5;

int date2mjd_nf(const char *date, double *mjd) {

  int nr, y, m, d, h, mm, s;
  long a,b;

  nr = sscanf( date, "%4d-%2d-%2dT%2d:%2d:%2d", &y, &m, &d, &h, &mm, &s);

  if ( nr == 6 ) {
    if ( m<=2 ) {
      y--;
      m += 12;
    }
    a = y/100;
    b = 2 - a + a/4;
    *mjd = b + floor(365.25 * y) + floor(30.6001 * (m + 1)) + d + 1720994.5;
//printf("mjd 1: %lf\n", *mjd);
//printf("%f\n", roundf((s + 60 * (mm + 60 * h)) / 86400. * 100000)/100000);
    *mjd += roundf((s + 60 * (mm + 60 * h)) / 86400. * 100000) / 100000;

    *mjd -= JD_OFF;

    return 0;

  } else {

    *mjd = 0.;
    return -1;

  }

  return 0;
}
