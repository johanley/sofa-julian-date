# sofa-julian-date

<a href='https://www.iausofa.org/'>SOFA</a> is a library from the IAU (International Astronomical Union) for fundamental astronomical calculations.

The SOFA library has functions for doing Julian date calculations.
Those functions are restricted to dates that are not before -4799 January 1 (in the Gregorian calendar).
Such restrictions are quite common for software that computes Julian date.

**This small project is an experiment in removing that restriction on dates.**

The algorithm used by SOFA for Julian date calculation is defined in the *Explanatory Supplement* (1992), <a href='https://archive.org/details/explanatorysuppl0003unse/page/604/mode/2up'>page 604</a>. The *Explanatory Supplement* (2006), <a href='https://archive.org/details/explanatorysuppl00pken/page/604/mode/2up'>page 604</a> has the same algorithm.
(Curiously, in the SOFA code, `iauCal2jd` is not superficially the same as the algo stated in the *Explanatory Supplement* 1992/2006.)

Fundamental astronomy has made giant strides in the past decades:
- measurements of stellar positions, proper motions, and radial velocities from instruments such as GAIA are remarkably precise.
- modern models of the Earth's precession and nutation provide accurate results stretching back tens of thousands of years.

Perhaps it's time to update the Julian date algorithm, and drop the restriction on dates.

## What I've Done

The SOFA code is completely intact, except for one small change: the `main` function in `t_sofa_c.c` has been renamed to `main_disabled`. 
This is simply because I want to implement my own `main` function for running unit tests.

The SOFA project uses C89, but I have used C99 in this project.

The files I have added are described below.

`alternate-jd-algos.c` defines two alternate implementations for two corresponding SOFA functions:
- alternate_iauCal2jd
- alternate_iauJd2cal

`alternate-headers.h` :
- standard C header file, with prototypes for the two alternate functions.

`alternate-run-tests.c` :
- defines and runs unit tests
- tests are run versus both the two original SOFA functions, and the two alternate implementations.

`run-tests.exe`
- a Windows executable that runs the tests ./run-tests.exe

`run-tests-results.txt`
- the output to `stdout` of running the tests.

## Bug Reports 

Bug reports about negative Julian dates in general:
-  https://github.com/astropy/astropy/issues/9231   
-  https://github.com/Unidata/netcdf4-python/issues/584
-  https://astronomy.stackexchange.com/questions/49790/calculation-of-julian-day-is-off-for-negative-dates

## Archeology

A 1968 paper by Fliegel and van Flandern in the *Communications of the ACM* is similar in nature to the algorithm stated in the *Explanatory Supplement* 1992/2006, but it's not the same (at least superficially):

https://dl.acm.org/doi/pdf/10.1145/364096.364097

