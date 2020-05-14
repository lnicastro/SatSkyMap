/*
mjd2date
   From MJD to ISO 8601 string date format (DATE-OBS), e.g. "2006-08-07T12:15:11".
   or
mjd2datef
   From MJD to date in the form "2006-08-07T12:15:11.123".


  Return '\0' for bad input date.

  LN@IASF-INAF, Aug 2006                       Last change: 30/01/2020
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* No fractional seconds */
char *mjd2date(double mjd)
{ 
  static char date[20];
  const double DJMIN = -68569.5;
  const double DJMAX = 1e9;
  const double dj1=2400000.5;

  int iy, im, id, h, mm, s;
  long jd, l, n, i, k;
  double dj, d1, d2, f1, f2, f, d;


/* Verify date is acceptable. */
  dj = dj1 + mjd;
  if (dj < DJMIN || dj > DJMAX) return 0;

/* Copy the date, big then small, and re-align to midnight. */
  if (dj1 >= mjd) {
     d1 = dj1;
     d2 = mjd;
  } else {
     d1 = mjd;
     d2 = dj1;
  }
  d2 -= 0.5;

/* Separate day and fraction. */
  f1 = fmod(d1, 1.);
  f2 = fmod(d2, 1.);
  f = fmod(f1 + f2, 1.);
  if (f < 0.) f += 1.;
  d = floor(d1 - f1) + floor(d2 - f2) + floor(f1 + f2 - f);
  jd = (long) floor(d) + 1L;

/* Express day in Gregorian calendar. */
  l = jd + 68569L;
  n = (4L * l) / 146097L;
  l -= (146097L * n + 3L) / 4L;
  i = (4000L * (l + 1L)) / 1461001L;
  l -= (1461L * i) / 4L - 31L;
  k = (80L * l) / 2447L;
  id = (int) (l - (2447L * k) / 80L);
  l = k / 11L;
  im = (int) (k + 2L - 12L * l);
  iy = (int) (100L * (n - 49L) + i + l);
  h = (int) (f*24);
  f -= h/24.;
  mm = (int) (f*1440);
  f -= mm/1440.;
  s = (int)(f*86400 + 0.5);

  sprintf(date, "%4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d", iy, im, id, h, mm, s);

  return date;
}


/* With fractional seconds (millisec) */
char *mjd2datef(double mjd)
{ 
  static char date[24];
  const double DJMIN = -68569.5;
  const double DJMAX = 1e9;
  const double dj1=2400000.5;

  int iy, im, id, h, mm, s, fs;
  long jd, l, n, i, k;
  double dj, d1, d2, f1, f2, f, d;


/* Verify date is acceptable. */
  dj = dj1 + mjd;
  if (dj < DJMIN || dj > DJMAX) return 0;

/* Copy the date, big then small, and re-align to midnight. */
  if (dj1 >= mjd) {
     d1 = dj1;
     d2 = mjd;
  } else {
     d1 = mjd;
     d2 = dj1;
  }
  d2 -= 0.5;

/* Separate day and fraction. */
  f1 = fmod(d1, 1.);
  f2 = fmod(d2, 1.);
  f = fmod(f1 + f2, 1.);
  if (f < 0.) f += 1.;
  d = floor(d1 - f1) + floor(d2 - f2) + floor(f1 + f2 - f);
  jd = (long) floor(d) + 1L;

/* Express day in Gregorian calendar. */
  l = jd + 68569L;
  n = (4L * l) / 146097L;
  l -= (146097L * n + 3L) / 4L;
  i = (4000L * (l + 1L)) / 1461001L;
  l -= (1461L * i) / 4L - 31L;
  k = (80L * l) / 2447L;
  id = (int) (l - (2447L * k) / 80L);
  l = k / 11L;
  im = (int) (k + 2L - 12L * l);
  iy = (int) (100L * (n - 49L) + i + l);
  h = (int) (f*24);
  f -= h/24.;
  mm = (int) (f*1440);
  f -= mm/1440.;
  s = (int)(f*86400);
  f -= s/86400.;
  fs = (int)(f*8.64e7 + 0.5);

  sprintf(date, "%4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%3.3d", iy, im, id, h, mm, s, fs);

  return date;
}
