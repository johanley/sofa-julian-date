#include "sofa.h"
#include "sofam.h"
#include <float.h>
#include <stdio.h>

/*
 Alternate algorithms for SOFA's iaucal2jd and iaujd2Cal functions, implemented in C99.

 These alternate algorithms aren't restricted to JD >= 0.

 Removing the JD >= 0 restriction now seems appropriate:
   - modern models for precession and nutation give remarkably precise results going back tens of thousands of years.
   - GAIA and other instruments yield very precise positions, proper motions, and radial velocities. 
   - Combined, these two allow the computation of stellar positions over extended time scales.

 Roughly speaking, these alternate implementations simply 'count the cycles', starting with largest cycles first.

 Both SOFA's implementations and these alternate implementations are for the Gregorian calendar.
*/

static const int NORMAL_YEAR_LEN = 365;
static const int LEAP_YEAR_LEN = 366;
static const int MONTH_LEN_TABLE[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const double JAN_0_YEAR_0 = 1721058.5; //January 0.0, of year 0 
static const int FULL_CYCLE_YEARS = 400;
static const int FULL_CYCLE_DAYS = (3*NORMAL_YEAR_LEN + LEAP_YEAR_LEN) * 25 * 4 /*centuries*/ - 3 /*oddball century-years with no leap day*/;

/* Leap year logic for the Gregorian calendar. */
static int is_leap(int y){
   int is_century_year = (y % 100 == 0);
   return is_century_year ? (y % 400 == 0) : (y % 4 == 0);
}

/* The number of days in the given year. */
static int year_len(int y) {
   return is_leap(y) ? LEAP_YEAR_LEN : NORMAL_YEAR_LEN;
}

/*
 Number of days in a full set of complete years.
 Includes the start-year, but excludes the end-year.
 Returns 0 if the start and end are the same year. 
*/
static int days_in_complete_years(int start_year, int end_year){
   int result = 0;
   for(int y = start_year; y < end_year; ++y) {
     result += year_len(y);
   }
   return result;
}

/* The length of the given month in days. The month-index is 1-based. */
static int month_len(int y, int m){
   int length = MONTH_LEN_TABLE[m - 1];
   int has_leap_day = is_leap(y) && (m == 2);
   return has_leap_day ? length + 1 : length;
}

/* 
 For the given date, return the number of days since Jan 0.0.
 Jan 0.0 is just an alias for December 31 of the previous year.
*/
static int days_from_jan0(int y, int m, int d){
  int days_in_completed_months = 0;
  //start with January, and move forward
  for(int completed_month = 1; completed_month < m; ++completed_month) {
    days_in_completed_months += month_len(y, completed_month); 
  }
  return days_in_completed_months + d;
}

/* The number of days remaining in the given month, from the given day. */
static int days_remaining_in_month(int y, int m, int d) {
  int length = month_len(y, m);
  return length + 1 - d;
}

/* Return the number of days until Dec 32.0 in the Gregorian calendar, for the given year and month. */
static int days_from_dec32(int y, int m, int d) {
   int days_in_completed_months = 0;
   //start with December, and count backwards in time
   for(int completed_month = 12; completed_month > m; --completed_month) {
     days_in_completed_months += month_len(y, completed_month); 
   }
   return days_in_completed_months + days_remaining_in_month(y, m, d);
}

static void cal2jd_non_neg_years(int iy, int im, int id, double *djm0, double *djm){
   // 1. full cycles in the Gregorian calendar 
   int num_cycles = iy / FULL_CYCLE_YEARS;
   int full_cycles = num_cycles * FULL_CYCLE_DAYS;

   // 2. remainder-years: whole years left after the full cycles
   int remainder_years = days_in_complete_years(num_cycles * FULL_CYCLE_YEARS, iy);
   
   // 3. remainder-days in the final year
   int remainder_days = days_from_jan0(iy, im, id);

   //REVISIT THIS LOGIC later
   double jd = JAN_0_YEAR_0 + full_cycles + remainder_years + remainder_days;
   *djm0 = jd;
   *djm = 0.0;
}

static void cal2jd_neg_years(int iy, int im, int id, double *djm0, double *djm) {
    //In the negative years, it's convenient to use (year + 1) as the base from which to track cycles.
    //This is because we're counting backwards through the calendar
    int y_biased = iy + 1;

    //1. full cycles in the calendar  
    int num_cycles = y_biased / FULL_CYCLE_YEARS;
    int full_cycles = abs(num_cycles * FULL_CYCLE_DAYS);
    
    //2. remainder years: whole years left after the full cycles
    int remainder_years = days_in_complete_years(y_biased, num_cycles * FULL_CYCLE_YEARS); 
    
    //3. remainder days in the final year
    int remainder_days = days_from_dec32(iy, im, id);
    int OVERHANG = 1; // Jan 0.0 is already impinging onto the negative years, by 1 day
    int total = full_cycles + remainder_years + remainder_days;

   //REVISIT THIS LOGIC later
   double jd = JAN_0_YEAR_0 + OVERHANG - total;
   *djm0 = jd;
   *djm = 0.0;
}

static void alternate_cal2jd(int iy, int im, int id, double *djm0, double *djm) {
   if (iy >= 0) {
      cal2jd_non_neg_years(iy, im, id, djm0, djm);
   }
   else {
      cal2jd_neg_years(iy, im, id, djm0, djm);
   }
}

// ***************************************
// AN ALTERNATE IMPLEMENTATION OF cal2jd.
// ***************************************
int alternate_iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
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
   if ( (id < 1) || (id > (MONTH_LEN_TABLE[im-1] + ly))) j = -3;

  alternate_cal2jd(iy, im, id, djm0, djm);
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



// ***************************************************
// BELOW IS FOR AN ALTERNATE IMPLEMENTATION OF jd2cal.
// ***************************************************

static void jd2cal_non_neg_years(double jd, int *iy, int *im, int *id, double *fd){
   double BASE = JAN_0_YEAR_0 + 1.0; 

    //1. full cycles in the calendar  
    double target = jd - BASE; //the target value we'll match below
    int num_full_cycles = (int)floor(target / FULL_CYCLE_DAYS);
    int year = num_full_cycles * FULL_CYCLE_YEARS; //starting value for the year; can increase below

    //this temp value is less than the target value, and approaches it from below
    double temp_target = num_full_cycles * FULL_CYCLE_DAYS; 

    //2. remainder years: whole years left after the full cycles (not including the final year)
    int year_full_cycles = year; //simply to remember this value in the loop below 
    for(int remainder_year_idx = 0; remainder_year_idx < FULL_CYCLE_YEARS; ++remainder_year_idx ) {
      double one_more_year = year_len(year_full_cycles + remainder_year_idx);
      if ((temp_target + one_more_year) <= target) {
        temp_target += one_more_year;
        ++year;
      } else { break; }
    }

    //3. months and days in the final year
    int month = 1; //January as starting point; can increase below
    for(int month_idx = 1; month_idx <= 12; ++month_idx) {
      double one_more_month = month_len(year, month_idx);
      if ((temp_target + one_more_month) <= target) {
        temp_target += one_more_month;
        ++month;
      } else { break; }
    }
    double days = target - temp_target + 1.0; //+1 since the base is Jan 1 0h, not Dec 31 0h

    *iy = year;
    *im = month;
    *id = dint(days); //rounds towards 0
    *fd = days - *id;
}

static void jd2cal_neg_years(double jd, int *iy, int *im, int *id, double *fd){
   double BASE = JAN_0_YEAR_0 + 1.0; 

   //1. full cycles in the calendar  1
   double target = jd - BASE; //the target value we'll match below
   int num_full_cycles = floor(target / FULL_CYCLE_DAYS) + 1 ;
   int year = num_full_cycles * FULL_CYCLE_YEARS; //starting value for the year; can decrease below
   --year; //because going backwards through the calendar

   //this temp value is more than the target value, and approaches it from above
   double temp_target = num_full_cycles * FULL_CYCLE_DAYS;

   //2. remainder years: whole years left after the full cycles (not including the final year)
   int year_full_cycles = year; //simply to remember this value in the loop below 
   for(int remainder_year_idx = 0; remainder_year_idx < FULL_CYCLE_YEARS; ++remainder_year_idx ) {
     double one_less_year = year_len(year_full_cycles - remainder_year_idx);
     if ((temp_target - one_less_year) > target) {
       temp_target -= one_less_year;
       --year;
     } else { break; }
   }

   //3. months and days in the final year
   int month = 12; //starting point; can decrease below
   for(int month_idx = 12; month_idx >= 1; --month_idx) { //go backwards, Dec to Jan!
     double one_less_month = month_len(year, month_idx);
     if ((temp_target - one_less_month) > target) {
       temp_target -= one_less_month;
       --month;
     } else { break; }
   }
    //count backwards from the end of the month
    double days = month_len(year, month) + 1 + target - temp_target;  //32 + (-0.5) = 31.5 for a time on Dec 31, for example 

    *iy = year;
    *im = month;
    *id = dint(days); //rounds towards 0
    *fd = days - *id;
}

static void alternate_jd2Cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd){
   double jan_1_year_0 = JAN_0_YEAR_0 + 1; 
   double jd = dj1 + dj2;
   if (jd >= jan_1_year_0) {
      jd2cal_non_neg_years(jd, iy, im, id, fd);
   }
   else {
      jd2cal_neg_years(jd, iy, im, id, fd);
   }
}

// ***************************************
// AN ALTERNATE IMPLEMENTATION OF jd2cal.
// ***************************************
int alternate_iauJd2cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd)
{
/* Minimum and maximum allowed JD */
   //const double DJMIN = -68569.5;  //NO LONGER NEEDED
   const double DJMAX = 1e9;

   long jd, i;
   double dj, f1, f2, d, s, cs, v[2], x, t, f;

/* Verify date is acceptable. */
   dj = dj1 + dj2;
   //if (dj < DJMIN || dj > DJMAX) return -1;
   if (dj > DJMAX) return -1;

/* Separate day and fraction (where -0.5 <= fraction < 0.5). */
   d = dnint(dj1);
   f1 = dj1 - d;
   jd = (long) d;
   d = dnint(dj2);
   f2 = dj2 - d;
   jd += (long) d;

/* Compute f1+f2+0.5 using compensated summation (Klein 2006). */
   s = 0.5;
   cs = 0.0;
   v[0] = f1;
   v[1] = f2;
   for ( i = 0; i < 2; i++ ) {
      x = v[i];
      t = s + x;
      cs += fabs(s) >= fabs(x) ? (s-t) + x : (x-t) + s;
      s = t;
      if ( s >= 1.0 ) {
         jd++;
         s -= 1.0;
      }
   }
   f = s + cs;
   cs = f - s;

/* Deal with negative f. */
   if ( f < 0.0 ) {

   /* Compensated summation: assume that |s| <= 1.0. */
      f = s + 1.0;
      cs += (1.0-f) + s;
      s  = f;
      f = s + cs;
      cs = f - s;
      jd--;
   }

/* Deal with f that is 1.0 or more (when rounded to double). */
   if ( (f-1.0) >= -DBL_EPSILON/4.0 ) {

   /* Compensated summation: assume that |s| <= 1.0. */
      t = s - 1.0;
      cs += (s-t) - 1.0;
      s = t;
      f = s + cs;
      if ( -DBL_EPSILON/2.0 < f ) {
         jd++;
         f = gmax(f, 0.0);
      }
   }

   //? Not sure if I should be doing anything with 'f' here.
   // At the moment, it is simply abandoned here.

   alternate_jd2Cal(dj1, dj2, iy, im, id, fd);
/* Express day in Gregorian calendar. */
   // l = jd + 68569L;
   // n = (4L * l) / 146097L;
   // l -= (146097L * n + 3L) / 4L;
   // i = (4000L * (l + 1L)) / 1461001L;
   // l -= (1461L * i) / 4L - 31L;
   // k = (80L * l) / 2447L;
   // *id = (int) (l - (2447L * k) / 80L);
   // l = k / 11L;
   // *im = (int) (k + 2L - 12L * l);
   // *iy = (int) (100L * (n - 49L) + i + l);
   // *fd = f;

/* Success. */
   return 0;

}