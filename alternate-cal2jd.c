#include "sofa.h"
#include "sofam.h"
#include <float.h>
#include <stdio.h>

/*
 Alternate algorithm for SOFA's iauCal2jd function, implemented in C99.

 This alternate algorithm doesn't fail for dates before -4799 January 1.

 Removing the JD date-restriction now seems appropriate:
   - modern models for precession and nutation give remarkably precise results going back tens of thousands of years.
   - GAIA and other instruments yield very precise positions, proper motions, and radial velocities. 
   - Combined, these two allow the computation of stellar positions over extended time scales.

 For the Gregorian calendar.
*/

static const double JAN_0_YEAR_0 = 1721058.5; //January 0.0 of year 0 = Dec 31.0 of year -1. 
static const int CYCLE_YEARS = 400;
static const int SHORT_YR = 365;
static const int LONG_YR = 366;
static const int MONTH_LEN[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

/*
 Convert a date in the Gregorian calendar to a Julian date.
 There's no restriction on the date.

 I base the calculation on counting days from January 0, year 0.
 Then I simply re-base the result at the end, to reflect the 
 usual origin-day for Julian dates.
 This exploits the (near) symmetry of the calendar cycles.

 I'm using a nice trick in Robin O'Leary's algorithm: 
   https://pdc.ro.nu/jd-code.html
*/
static void gregorian_cal_to_julian_date(int y, int m, double d, double *djm0, double *djm)  {
    //completed years: small asymmetry between positive and negative years
    int y_p = (y >= 0) ? (y - 1) : y;  //y_p = y-prime
    int num_366yrs = (y_p/4) - (y_p/100) + (y_p/CYCLE_YEARS); //Robin's clever trick
    if (y > 0) {
      num_366yrs += 1; //since year 0 is a leap year
    }
    int num_365yrs = y - num_366yrs;
    double res = num_365yrs * SHORT_YR + num_366yrs * LONG_YR;    
    
    //completed months
    //Explanatory Supplement 1961, page 434: 
    int DAYS_IN_PRECEDING_MONTHS[] = {0 /*Jan*/, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /*Dec*/};
    res += DAYS_IN_PRECEDING_MONTHS[m-1];   
    int is_leap = (y % 100 == 0) ? (y % 400 == 0) : (y % 4 == 0);
    res += (is_leap && (m - 1) >= 2 ? 1 : 0); //'correct' for leap years  
    
    res += d;  // the day of the month
    
    //rebase to the usual origin of Julian date
    res += JAN_0_YEAR_0;   

   //REVISIT THIS LOGIC later
   *djm0 = res;
   *djm = 0.0;
}



// AN ALTERNATE IMPLEMENTATION OF cal2jd, which calls the function above instead.
int terse_alternate_iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
{
   int j, ly;

   //NO LONGER NEEDED
/* Earliest year allowed (4800BC) */
   //const int IYMIN = -4799; 

   //THIS IS USED THIS IN TWO PLACES NOW; MOVED TO A STATIC, NEAR TOP OF FILE
/* Month lengths in days */
   //static const int mtab[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

/* Preset status. */
   j = 0;

/* Validate year and month. */
   //if (iy < IYMIN) return -1; //NO LONGER NEEDED
   if (im < 1 || im > 12) return -2;

/* If February in a leap year, 1, otherwise 0. */
   ly = ((im == 2) && !(iy%4) && (iy%100 || !(iy%400)));

/* Validate day, taking into account leap years. */
   if ( (id < 1) || (id > (MONTH_LEN[im-1] + ly))) j = -3;

  gregorian_cal_to_julian_date(iy, im, id, djm0, djm);
/* Return result. */
   // my = (im - 14) / 12;
   // iypmy = (long) (iy + my);
   // *djm0 = DJM0;
   // *djm = (double)((1461L * (iypmy + 4800L)) / 4L
   //               + (367L * (long) (im - 2 - 12 * my)) / 12L
   //               - (3L * ((iypmy + 4900L) / 100L)) / 4L
   //               + (long) id - 2432076L);

/* Return status. */
   return j;

/* Finished. */

}

