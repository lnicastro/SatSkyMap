/*
  MJD now from system time.

  LN@INAF-OAS, Nov. 2020                      Last change: 08/11/2020
*/

#include <stdio.h>
#include <strings.h>
#include <time.h>


// Julian epoch MJD
double getMJD_J(int day, int month, int year, int hour, int minute, int second) {

  int a = (14 - month) / 12;
  int y = year + 4800 - a;
  int m = month + 12*a - 3;
   
  double d = day + (hour + minute/60e0 + second/3600e0) / 24.;
 
  return d + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045 - 2400000;
}      


int main (void) {

  time_t now = time(&now);
  struct tm *ptm = gmtime(&now);

//char buf[256] = {0};
//strftime(buf, 256, "%d/%m/%Y", ptm);
//puts(buf);
    
  int year = ptm->tm_year + 1900;
  int month = ptm->tm_mon + 1;
  int day = ptm->tm_mday;
  int hour = ptm->tm_hour;
  int minute = ptm->tm_min;
  int second = ptm->tm_sec;

  double mjd = getMJD_J(day, month, year, hour, minute, second);

  printf("%lf\n", mjd); 

  return 0;
}
