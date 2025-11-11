#include <stdio.h>
#include <math.h>
#include "sofa.h"
#include "alternate-headers.h"

/*
 UNIT TESTS FOR iauCal2jd and iauJd2cal from SOFA, and also for alternate implementations to 
 those same two functions.

 Implemented in C99.

 In PowerShell on Windows, in the current directory:
   gcc $(Get-ChildItem -Path *.c -Name) -std=c99 -Wall -Werror -Wpedantic -o run_tests.exe
   ./run_tests.exe
*/

static int num_errors = 0;
static int num_successful = 0;

static const char *SUCCESS = "OK";
static const char *FAILURE = " X";

static const char *SOFA = "SOFA";
static const char *ALTERNATIVE = "ALT ";

static const int REPORT = 1;
static const int SILENT = 0;

/* Output the result of checking expected-result versus actual-result. */
static void check_date_to_jd(double expected, double result, const char *source, const int report){
    const char* message = (expected == result) ? SUCCESS : FAILURE;
    short ok = expected == result;
    !ok ? ++num_errors : ++num_successful;
    if (report){
      printf("%s %s Expected: %f Result: %f\n", source, message, expected, result);
    }
}

/* Output the result of checking expected-result versus actual-result. */
static void check_jd_to_date(int expected_y, int expected_m, int expected_d, double expected_fd, 
    int result_y, int result_m, int result_d, double result_fd, const char *source, const int report){
    short ok_fd = abs(expected_fd - result_fd) < __FLT_EPSILON__;
    short ok = (expected_y == result_y) &&  (expected_m == result_m) && (expected_d == result_d) && ok_fd;
    const char* message_a = ok ? SUCCESS : FAILURE;
    !ok ? ++num_errors : ++num_successful;
    if (report){
        printf("%s %s Expected: %d-%d-%d %f Result: %d-%d-%d %f\n", \
            source,
            message_a, 
            expected_y, expected_m, expected_d, expected_fd, \
            result_y, result_m, result_d, result_fd
        );
    }
}

static void report_error(int status){
    if (status != 0){
      printf("Error status: %d\n", status);
    }
}

/* Test converting a date in the Gregorian calendar to a Julian date. */
static void test_yyyy_mm_dd_to_julian_date(
    int y, int m, int d, double fd, double expected_jd, 
    int (*func)(int, int, int, double*, double*), const char *source, int report
){
    //there's no fractional-day input here; this is asymmetric with respect to the complementary function, 
    //which returns a fractional-day. This means I need to add a fractional day 'manually', after the calc for the whole day
    double djm0 = 0.0;
    double djm = 0.0;
    int status = func(y, m, d, &djm0, &djm); 
    report_error(status);
    double result = djm0 + djm + fd;
    check_date_to_jd(expected_jd, result, source, report);
}

/** Test converting a Julian date to a date in the Gregorian calendar. */
static void test_julian_date_to_yyyy_mm_dd(
    double jd1, double jd2, 
    int expected_y, int expected_m, int expected_d, double expected_fd, 
    int (*func)(double, double, int*, int*, int*, double*), const char *source, int report
){
    int result_y, result_m, result_d;
    double result_fd;
    int status = func(jd1, jd2, &result_y, &result_m, &result_d, &result_fd); 
    report_error(status);
    check_jd_to_date(expected_y, expected_m, expected_d, expected_fd, result_y, result_m, result_d, result_fd, source, report);
}

/* 
 Test the conversions in both directions, jd-to-calendar-date, and calendar-date-to-jd. 
 Test both the original SOFA function, and then the alternate implementation.
 Thus, a call to this function performs 4 tests in sequence.
*/ 
static void test_both_directions(int y, int m, int d, double fd, double jd1, double jd2){
    test_julian_date_to_yyyy_mm_dd(jd1, jd2, y, m, d, fd, iauJd2cal, SOFA, REPORT);
    test_julian_date_to_yyyy_mm_dd(jd1, jd2, y, m, d, fd, alternate_iauJd2cal, ALTERNATIVE, REPORT);
    test_yyyy_mm_dd_to_julian_date(y, m, d, fd, jd1 + jd2, iauCal2jd, SOFA, REPORT);
    test_yyyy_mm_dd_to_julian_date(y, m, d, fd, jd1 + jd2, alternate_iauCal2jd, ALTERNATIVE, REPORT);
    printf("\n");
}

/* These tests aren't reported in detail. Only the count of success-fail is reported for these. */ 
static void test_entire_year(int y, double jd_jan_0){
    printf("Testing every day in the year: %d\n", y);
    int day_num = 0; //1..(365|366)
    for(int m = 1; m <= 12; ++m){
        int num_days_in_month = month_len(y, m);
        for(int d = 1; d <= num_days_in_month; ++d){
            ++day_num;
            double jd = jd_jan_0 + day_num;
            test_julian_date_to_yyyy_mm_dd(jd, 0.0, y, m, d, 0.0, iauJd2cal, SOFA, SILENT);
            test_julian_date_to_yyyy_mm_dd(jd, 0.0, y, m, d, 0.0, alternate_iauJd2cal, ALTERNATIVE, SILENT);
            test_yyyy_mm_dd_to_julian_date(y, m, d, 0.0, jd, iauCal2jd, SOFA, SILENT);
            test_yyyy_mm_dd_to_julian_date(y, m, d, 0.0, jd, alternate_iauCal2jd, ALTERNATIVE, SILENT);
        }
    }
}

/* Test every day of the year for years near the year 0. (These cases are easy to calculate manually.) */ 
static void test_small_years() {
    double base = 1721058.5;
    test_entire_year(-9, base - 2*366 - 7*365);
    test_entire_year(-8, base - 2*366 - 6*365);
    test_entire_year(-7, base - 1*366 - 6*365);
    test_entire_year(-6, base - 1*366 - 5*365);
    test_entire_year(-5, base - 1*366 - 4*365);
    test_entire_year(-4, base - 1*366 - 3*365);
    test_entire_year(-3, base - 0*366 - 3*365);
    test_entire_year(-2, base - 0*366 - 2*365);
    test_entire_year(-1, base - 0*366 - 1*365);
    test_entire_year(0, base + 0*366 + 0*365);
    test_entire_year(1, base + 1*366 + 0*365);
    test_entire_year(2, base + 1*366 + 1*365);
    test_entire_year(3, base + 1*366 + 2*365);
    test_entire_year(4, base + 1*366 + 3*365);
    test_entire_year(5, base + 2*366 + 3*365);
    test_entire_year(6, base + 2*366 + 4*365);
    test_entire_year(7, base + 2*366 + 5*365);
    test_entire_year(8, base + 2*366 + 6*365);
    test_entire_year(9, base + 3*366 + 6*365);
    test_entire_year(10, base + 3*366 + 7*365);
    test_entire_year(11, base + 3*366 + 8*365);
    test_entire_year(12, base + 3*366 + 9*365);
  }


/* 
 Run all tests for conversions from calendar-date to Julian date, and vice versa. 
 This tests both the original SOFA algorithms, and the 2 alternate algorithms implemented in alternate-algos.c.
*/
static void run_tests_for_both_old_and_new_algorithms(){
    printf("SOFA's tests.\n");
    test_both_directions(2003, 6, 1, 0.0, 2400000.5, 52791.0);

    //test_julian_date_to_yyyy_mm_dd(2400000.5, 50123.9999, 1996, 2, 10, 0.9999);

    test_both_directions(1996, 2, 11, 0.0, 2400000.5, 50124.0);  //my modification of SOFA's test, in order to use whole days
    /*
    FROM t_sofa_c.c:

    dj1 = 2400000.5;
    dj2 = 50123.9999;
    j = iauJd2cal(dj1, dj2, &iy, &im, &id, &fd);
    viv(iy, 1996, "iauJd2cal", "y", status);
    viv(im, 2, "iauJd2cal", "m", status);
    viv(id, 10, "iauJd2cal", "d", status);
    vvd(fd, 0.9999, 1e-7, "iauJd2cal", "fd", status);
    viv(j, 0, "iauJd2cal", "j", status);

    j = iauCal2jd(2003, 06, 01, &djm0, &djm);
    vvd(djm0, 2400000.5, 0.0, "iauCal2jd", "djm0", status);
    vvd(djm,    52791.0, 0.0, "iauCal2jd", "djm", status);
    */

    printf("\nExplanatory Supplement, 1961, page  437.\n");
    test_both_directions(1500, 1, 1, 0.0, 2268923.5, 0.0);
    test_both_directions(1600, 1, 1, 0.0, 2305447.5, 0.0);
    test_both_directions(1700, 1, 1, 0.0, 2341972.5, 0.0);
    test_both_directions(1800, 1, 1, 0.0, 2378496.5, 0.0);
    test_both_directions(1900, 1, 1, 0.0, 2415020.5, 0.0);

    test_both_directions(1500, 3, 1, 0.0, 2268923 + 0.5 + 59, 0.0); 
    test_both_directions(1600, 3, 1, 0.0, 2305447 + 0.5 + 60, 0.0); //March 1 is after Feb 29; only this one is a leap year
    test_both_directions(1700, 3, 1, 0.0, 2341972 + 0.5 + 59, 0.0);
    test_both_directions(1800, 3, 1, 0.0, 2378496 + 0.5 + 59, 0.0);
    test_both_directions(1900, 3, 1, 0.0, 2415020 + 0.5 + 59, 0.0);

    printf("\nGuide de Donnees Astronomiques 2017, Bureau des longitudes, page 8.\n");
    test_both_directions(1950, 1, 1, 0.5, 2433283.0, 0.0);
    test_both_directions(2000, 1, 1, 0.5, 2451545.0, 0.0);
    test_both_directions(2050, 1, 1, 0.5, 2469808.0, 0.0);
    test_both_directions(2090, 1, 1, 0.5, 2484418.0, 0.0);

    // -1374 May 3, at 13:52:19.2 TT 
    printf("\nFrom Vondrak, Wallace, Capitaine 2011.\n");
    test_both_directions(-1374, 5, 3, 0.578, 1219339.078, 0.0);    

    printf("\nObserver's Handbook, RASC, 2024, page 47.\n");
    test_both_directions(2024, 1, 1, 0.0, 2460310.5, 0.0);    
    test_both_directions(2024, 3, 1, 0.0, 2460370.5, 0.0);    

    printf("\nAstronomical Algorithms, Meeus 1991, page 61ff.\n");
    test_both_directions(1957, 10, 4, 0.81, 2436116.31, 0.0);
    test_both_directions(1987, 6, 19, 0.5, 2446966.0, 0.0);

    printf("\nFrom https://legacy-www.math.harvard.edu/computing/javascript/Calendar/index.html\n");
    test_both_directions(-8, 1, 1, 0.5, 1718138.0, 0.0);
    test_both_directions(-101, 1, 1, 0.5, 1684171.0, 0.0);
    test_both_directions(-799, 1, 1, 0.5, 1429232.0, 0.0);
    test_both_directions(-800, 1, 1, 0.5, 1428866.0, 0.0);
    test_both_directions(-801, 1, 1, 0.5, 1428501.0, 0.0);
    test_both_directions(99, 12, 31, 0.5, 1757584.0, 0.0);
    test_both_directions(100, 1 , 1, 0.5, 1757584.0 + 1.0, 0.0);
    test_both_directions(100, 1, 31, 0.5, 1757584.0 + 31.0, 0.0);
    test_both_directions(100, 2, 1, 0.5, 1757584.0 + 31.0 + 1.0, 0.0);
    test_both_directions(100, 2, 28, 0.5, 1757584.0 + 31.0 + 28.0, 0.0); //100 is not a leap year
    test_both_directions(100, 3, 1, 0.5, 1757584.0 + 31.0 + 28.0 + 1.0, 0.0);
    test_both_directions(3000, 1, 1, 0.5, 2816788, 0.0);
    test_both_directions(30000, 1, 1, 0.5, 12678335.0, 0.0);
    test_both_directions(100, 1, 1,0.5, 1757585.0, 0.0);
    test_both_directions(101, 1, 1, 0.5, 1757950.0, 0.0); 
    test_both_directions(200, 1, 1, 0.5, 1794109.0, 0.0); 
    test_both_directions(300, 1, 1, 0.5, 1830633.0, 0.0); 
    test_both_directions(400, 1, 1, 0.5, 1867157, 0.0); 
    test_both_directions(700, 1, 1, 0.5, 1976730, 0.0);  
    test_both_directions(800, 1, 1, 0.5, 2013254, 0.0);

    printf("\nThe origin of the Julian date is -4712-01-01 12h, in the Julian calendar. That is -4713-11-24 in the Gregorian calendar.\n");
    test_both_directions(-4713, 11, 24, 0.5, 0.0, 0.0);

    printf("\nThe first date supported by the SOFA algorithm: -4799-01-01.\n");
    test_both_directions(-4799, 1, 1, 0.0, -31738.5, 0.0);

    printf("\nNum failed tests: %d\n", num_errors);
    printf("Num successful tests: %d\n", num_successful);

    printf("\nTest entire years near the year 0.\n");
    printf("There's no detailed reporting in these cases.\n");
    test_small_years();

    printf("\nNum failed tests: %d\n", num_errors);
    printf("Num successful tests: %d\n", num_successful);
}

/* (I have renamed the 'main' function found in t_sofa_c.c, in order to replace it with this 'main'.)  */
int main(){
  run_tests_for_both_old_and_new_algorithms();
}