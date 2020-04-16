/*

  Define useful constants

*/

#ifndef CONST_DEF_H
#define CONST_DEF_H

/* From MJD to JD */
#define JD_OFF 2400000.5

/* Julian date at standard epoch 1.5 Jan 2000 */
#define JD2000 2451545.0

/* degrees to radians */
#define DEG2RAD 0.017453292519943295769236907684886127134428718885417
//#define DD2R 1.74532925199432957692369E-2

/* radians to degrees */
#define RAD2DEG 57.295779513082320876798155

/* hours to radians */
#define HRS2RAD 0.26179938779914943653855361527329190701643078328126

/* radians to hours */
#define RAD2HRS 3.8197186342054880584532103209403446888270314977709


/* PI and 2PI */
#define PI 3.141592653589793238462643383279
#define TWOPI 6.2831853071795864769252867665590057683943387987502

/* seconds in a day */
#define SEC_IN_DAY 86400.

#define TIME_EPSILON (1./SEC_IN_DAY)  // 1s

#endif
