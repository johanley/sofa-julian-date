## sofa-julian-date

<a href='https://www.iausofa.org/'>SOFA</a> is a library from the IAU (International Astronomical Union) for fundamental astronomical calculations.

The SOFA library has functions for doing Julian date calculations.
Those functions are restricted to dates that are not before -4799 January 1 (in the Gregorian calendar).
Such restrictions are quite common for software that computes Julian date.

This small project is an experiment in removing that restriction. 

The algorithm used by SOFA for Julian date calculation is defined in the *Explanatory Supplement*, 2nd edition, 2006, <a href='https://archive.org/details/explanatorysuppl00pken/page/606/mode/2up'>page 606</a>.
That algorithm dates from 1970.

Fundamental astronomy has made giant strides since 1970:
- measurements of stellar positions, proper motions, and radial velocities from instruments such as GAIA are remarkably precise.
- modern models of the Earth's precession and nutation provide accurate results stretching back tens of thousands of years.

Perhaps it's time to update the Julian date algorithm, and drop the restriction on dates.

# What I've Done

The SOFA code is completely intact, except for one small change: the `main` function in `t_sofa_c.c` has been renamed to `main_disabled`. 
This because I want to implement my own `main` function for running my own unit tests.

The SOFA project uses C89. I have used C99 in this project.

The files I have added are described below.

`alternate-jd-algos.c` defines two alternate implementations for two corresponding SOFA functions:
- alternate_iauCal2jd
- alternate_iauJd2cal

`alternate-headers.h` :
- standard C header file, with prototypes for the two alternate functions.

`alternate-run-tests.c` :
- unit tests
- tests are run versus both the two original SOFA functions, and the new alternate implementations.

run-tests.exe
- a Windows executable that runs the tests ./run-tests.exe

run-tests-results.txt
- the output to `stdout` of running the tests.

