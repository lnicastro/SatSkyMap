/*
   From DATE in the form "2006-08-07 12:15:11" or "2006-08-07 12:15:11.123" to MJD.
   (from 1 to 3 fractional digits are allowed)

   Note the " " (ISO 8601 date format) can be "T" (FITS).
   Missing "time" info is managed too.

  Return 0 un success.

  LN@INAF-OAS, Aug 2006                       Last change: 30/01/2020
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int date2mjd(const char *date, double *mjd) {

  char *chr;
  double fs_fac=1e3;


  if (strlen(date) < 10)
    return -1;

/* Copy into a local buffer */
  char *d_o;
  d_o = (char *) calloc(strlen(date)+1, sizeof(char));
  strcpy(d_o,date);

// Eventually replace the T with blanck
  chr = strchr(d_o, 'T');
  if (chr != NULL) *chr = ' ';

  *mjd = 0.;

/* Split DATE-OBS (first element in STR kwds) into its components for the
   table columns */
  int y=0, m, d, h=0, mm=0, s=0, fs=0;

  chr = strchr(d_o, '.');
  if (chr == NULL || strlen(d_o) == 20)  // no fractional sec or just "." 
    sscanf(d_o, "%4d-%2d-%2d %2d:%2d:%2d", &y, &m, &d, &h, &mm, &s);
  else {
    if (strlen(d_o) >= 23)  // 3 fractional digits
      sscanf(d_o, "%4d-%2d-%2d %2d:%2d:%2d.%3d", &y, &m, &d, &h, &mm, &s, &fs);
    else if (strlen(d_o) == 22) {  // 2 fractional digits
      sscanf(d_o, "%4d-%2d-%2d %2d:%2d:%2d.%2d", &y, &m, &d, &h, &mm, &s, &fs);
      fs_fac = 1e2;
    } else if (strlen(d_o) == 21) {  // 1 fractional digits
      sscanf(d_o, "%4d-%2d-%2d %2d:%2d:%2d.%1d", &y, &m, &d, &h, &mm, &s, &fs);
      fs_fac = 1e1;
    } else if (strlen(d_o) == 16) {  // no secs
      sscanf(d_o, "%4d-%2d-%2d %2d:%2d", &y, &m, &d, &h, &mm);
      fs_fac = 1e1;
    } else if (strlen(d_o) == 13) {  // no mins and secs
      sscanf(d_o, "%4d-%2d-%2d %2d", &y, &m, &d, &h);
      fs_fac = 1e1;
    } else if (strlen(d_o) == 10) {  // no hour, mins and secs
      sscanf(d_o, "%4d-%2d-%2d", &y, &m, &d);
      fs_fac = 1e1;
    }
  }


  if (y < 0)
    y++;
  else if (y == 0)
    return 1;

  double day = d + ( h + mm/60.e0 + (s+fs/fs_fac)/3600.e0) / 24.;

  if ( m < 3 ) {  // Jan or Feb -> no leap day
    y--;
    m += 12;
 }

 long a = (long)y/100;

 *mjd = (long)(y*0.25e0) + 365.0e0*(y -1860.e0) + (long)(30.6001e0*(m+1.)) +
      day  - 106e0;

// Gregorian Calendar starts on Oct. 15, 1582 (MJD -100830.)
 if (*mjd > -100830.)
   *mjd += 2 - a + (long)(a/4);

  return 0;
}
