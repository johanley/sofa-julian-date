#include "sofa.h"
#include "sofam.h"
#include <float.h>
#include <stdio.h>

/*
 Alternate algorithm for SOFA's iauJd2cal function, implemented in C99.

 This alternate algorithm doesn't fail for dates before -4799 January 1.

 Removing the JD date-restriction now seems appropriate:
   - modern models for precession and nutation give remarkably precise results going back tens of thousands of years.
   - GAIA and other instruments yield very precise positions, proper motions, and radial velocities. 
   - Combined, these two allow the computation of stellar positions over extended time scales.

 For the Gregorian calendar.
*/

static const int SHORT_YR = 365;
static const int LONG_YR = 366;
static const double JAN_1_YEAR_0 = 1721058.5 + 1.0;
static const int CYCLE_YEARS = 400;
static const int CYCLE_DAYS = SHORT_YR*CYCLE_YEARS + CYCLE_YEARS/4 - CYCLE_YEARS/100 + CYCLE_YEARS/CYCLE_YEARS; //146_097 days
static const int MONTH_LEN[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

static int is_leap(int y) {
  return (y % 100 == 0) ? (y % 400 == 0) : (y % 4 == 0);
}

/* 
 The length of the given month in days. The month-index is 1-based. 
 This function is non-static, since it's used in the testing file.
*/
int the_month_len(int y, int m){
   int length = MONTH_LEN[m - 1];
   int has_leap_day = is_leap(y) && (m == 2);
   return has_leap_day ? length + 1 : length;
}

/*
  There's no restriction on the date.

  Mental model: use a 'base', a point in time occurring once every 400 years, 
  at which the calendar cycle starts.  Counting forward in time from such any 
  such 'base' exploits the symmetry of the calendar's cycle.

  Take a base as always falling on a N*400 years from January 1.0, year 0:
    JD of a base = 1_721_059.5 + N * 146_097, with  N = ...-2,-1,0,1,2,...

  There are 2 loops in this implementation, with a max number of 14 loop iterations.
*/
static void julian_date_to_gregorian_cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd) {
    double jd = dj1 + dj2;
    //1. find the closest 'base' that PRECEDES the given moment
    int num_cycles = (int)floor((jd - JAN_1_YEAR_0)/CYCLE_DAYS); //rounds towards negative infinity: good!
    double base_jd = JAN_1_YEAR_0 + num_cycles * CYCLE_DAYS; //a January 1.0 in the years  .., -400, 0, 400, ..
    int year = num_cycles * CYCLE_YEARS; // ..,-400, 0, 400,.. (the starting value)
    double jd_minus_base = jd - base_jd; //never neg

    //THE GAME IS: to move this cursor forward from our base Jan 1.0 to the target jd
    double cursor = 0.0; 
    
    //2. remainder-years: whole, completed years after the base 
    //one big chunk of years: calculate a MINIMUM number of full remainder-years, to reduce loop iterations below
    int approx_days = (int)floor(jd_minus_base);
    int more_years = (approx_days / LONG_YR) - 1; // at least this many
    if (more_years > 0) {
      int m_p = more_years - 1;
      int more_days = more_years * SHORT_YR + (m_p/4) - (m_p/100) + (m_p/400) + 1;
      cursor += more_days; //still on a Jan 1.0!
      year += more_years;
    }
    //loop to find the rest of the remaining-years: at most 2 iterations here!
    int year_so_far = year; //for use in the loop 
    for(int more = 0; more < CYCLE_YEARS; ++more ) { 
      int year_length = is_leap(year_so_far + more) ? LONG_YR : SHORT_YR;
      if (cursor + year_length <= jd_minus_base) {
        cursor += year_length; // Jan 1.0 of the next year
        ++year;
      } else { break; }
    }
    
    //3. months and days
    int month = 0; //both a loop index AND a result-value
    double fractional_days = 0.0;
    for(month = 1; month <= 12; ++month) {
      int month_length = MONTH_LEN[month - 1];
      if (is_leap(year) && month == 2) ++month_length;
      if (cursor + month_length <= jd_minus_base) {
        cursor += month_length; //1st day of the next month
      }
      else {
        fractional_days = jd_minus_base - cursor + 1.0; break;
      }
    }
    *iy = year;
    *im = month;
    double integer_part = 0.0;
    *fd = modf(fractional_days, &integer_part);
    *id = (int)integer_part;
}



// AN ALTERNATE IMPLEMENTATION OF jd2cal, which calls the function above instead.
int terse_alternate_iauJd2cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd)
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

   julian_date_to_gregorian_cal(dj1, dj2, iy, im, id, fd);
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