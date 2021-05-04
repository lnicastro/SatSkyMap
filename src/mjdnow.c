/*
  MJD now from system time.

  LN@INAF-OAS, Nov. 2020                      Last change: 04/05/2020
*/

#include <stdio.h>
#include <strings.h>
#include <time.h>


// Julian epoch MJD - post Oct. 15, 1582 (MJD -100830.)
double getMJD_J(int day, int month, int year, int hour, int minute, int second) {

  //int a = (14 - month) / 12;
  //int y = year + 4800 - a;
  //int m = month + 12*a - 3;
   
  if ( month < 3 ) {  // Jan or Feb -> no leap day
	year--;
	month += 12;
  }

  double d = day + (hour + minute/60e0 + second/3600e0) / 24.;
 
  return d + (long)(year*0.25e0) + 365e0*(year -1860e0) + (long)(30.6001e0*(month+1.)) - 106e0 + 2 - year/100 + year/400;
}      


int main (void) {

  time_t now = time(&now);
//struct tm *ptm = gmtime(&now);

//char buf[256] = {0};
//strftime(buf, 256, "%d/%m/%Y %H:%M:%S", ptm);
//puts(buf);
    
//int year = ptm->tm_year + 1900;
//int month = ptm->tm_mon + 1;
//int day = ptm->tm_mday;
//int hour = ptm->tm_hour;
//int minute = ptm->tm_min;
//int second = ptm->tm_sec;
//printf("%d %d %d %d %d %d\n", day, month, year, hour, minute, second);
//double mjd = getMJD_J(day, month, year, hour, minute, second);
//double mjd = getMJD_J(ptm->tm_mday, ptm->tm_mon + 1, ptm->tm_year + 1900, ptm->tm_hour, ptm->tm_min, ptm->tm_sec);

//printf("%lf\n", mjd);

printf("%lf\n", (double)(40587. + (double)(now / (double) 86400.)));

  return 0;
}
