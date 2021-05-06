/*
  MJD now from system time.

  LN@INAF-OAS, Nov. 2020                      Last change: 04/05/2020
*/

#include <stdio.h>
#include <strings.h>
#include <time.h>

int main (void) {

  time_t now = time(&now);
  printf("%lf\n", (double)(40587. + (double)(now / (double) 86400.)));

  return 0;
}
