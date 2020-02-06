/* Implementation of solar_geometry.
   Copyright (C) 1997-2020 MINES ParisTech
   This file is part of the solar_geometry library.

   The solar_geometry library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   The solar_geometry library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the solar_geometry library; if not, see
   <http://www.gnu.org/licenses/>.  */

/*-------------------------------------------------------------------------*/
/*                        ECOLE DES MINES DE PARIS                         */
/*        CENTRE D'ENERGETIQUE - GROUPE TELEDETECTION & MODELISATION       */
/*                       Rue Claude Daunesse, BP 207                       */
/*                   06904 Sophia Antipolis cedex, FRANCE                  */
/*          Tel (+33) 04 93 95 74 49     Fax (+33) 04 93 95 75 35          */
/*                       E-mail : (name)@cenerg.cma.fr                     */
/*-------------------------------------------------------------------------*/
/*   L. Wald - O. Bauer - February 1997                                    */
/*   modified 8 July 2004 L. Wald for geocentric - geographic lat          */
/*   modified Nov. 2011 by P. Blanc : Add POSITION OF THE SUN IN THE SKY   */
/*            FAST                                                         */
/*   modified 2019 by B. Gschwind : API cleanup & fixes                    */
/*-------------------------------------------------------------------------*/

#define __C_solar_geometry

#include "solar_geometry.h"

#include <math.h>
#include <stdio.h>

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))

#define SG1_PI_2     1.57079632679489661923

/*************/
/* NOTATIONS */
/*************/
/* phi_g  : geographic latitude of the site, positive to North */
/* phi    : geocentric latitude of the site, positive to North */
/* lambda : longitude of the site, positive to East */
/* delta  : solar declination angle */
/* omega  : solar hour angle */
/* gamma  : solar altitude (or elevation) angle */
/* theta  : solar incidence (or zenithal) angle (ESRA --> zeta) */
/* alpha  : solar azimuthal angle (or psi) */
/* t   : solar time = true solar time (TST) = local apparent time (LAT) */
/* LAT : local apparent time or solar time or true solar time (TST)
 --> this system of time offers the advantage of symmetry of the solar
 geometry about the north-south line */
/* LMT : local mean time or clock time */
/* UT  : Universal Time, is GMT measured from Greenwich mean midnight */
/* omega_sr : sunrise hour angle */
/* omega_ss : sunset hour angle */
/* t_sr     : time of astronomical sunrise */
/* t_ss     : time of astronomical sunset */
/* omega1 : solar hour angle at beginning of the time period */
/* omega2 : solar hour angle at end of the time period */

/* S0  : astronomical daylength or astronomical sunshine duration */
/* I0  : solar constant = annual mean value of extraterrestrial direct solar
 irradiance G0 (1367.0 W/m2) */
/* G0  : extraterrestrial global solar irradiation (on an horizontal plane)
 =B0 */
/* G0h : hourly extraterrestrial solar irradiation (on an horizontal plane) */
/* G0d : daily extraterrestrial solar irradiation (on an horizontal plane) */

/* NB : All angles are computed in radians as a standard.
 The basic trigonometric calculations on the position of the sun are
 carried out in LAT. */

/*********************************************/
/*                                           */
/* G E O M E T R Y  O F  S O L A R  B E A M, */
/*                                           */
/* A  N E A R  P O I N T  S O U R C E        */
/*                                           */
/*********************************************/

/********************/
/*
 * BASIC PARAMETERS 
 */
/********************/

/* the last value is used for sg1_day_of_year_to_ymd */
static int const SG1_MONTHLY_DAY_OF_YEAR_OFFSET[13] = {
        0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
};

/* the last value is used for sg1_day_of_year_to_ymd */
static int const SG1_MONTHLY_DAY_OF_YEAR_OFFSET_LEAP_YEAR[13] = {
        0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366
};

static int SG1_DAYS_PER_MONTH[12] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

inline static int is_leap_year(int const year) {
    return (((year % 4) == 0) && ((year % 100) != 0)) || ((year % 400) == 0);
}

/*
 * Source : 
 */
/*
 * Inputs : day_of_month : day of the month (1..31) month_number : month number (1..12)
 * year_number : year number (4 digits) 
 */
/*
 * Outputs : julian_day : integer day number or julian day (1..366) 
 */
/*
 * The procedure "make_julian_day" converts a day given in day, month and year into a
 * julian day. Returns 0 if OK, 1 otherwise. 
 */
int sg1_ymd_to_day_of_year(int year, int month, int day_of_month,
        int * const day_of_year)
{
    int julien;

    if (year <= 0)
        return 1;
    if ((month < 1) || (month > 12))
        return 1;

    /* technicaly this is not required */
    if ((day_of_month < 1) || (day_of_month > 31))
        return 1;

    *day_of_year = day_of_month + SG1_MONTHLY_DAY_OF_YEAR_OFFSET[month - 1];
    if (is_leap_year(year) && (month > 2))
        *day_of_year += 1;
    return 0;
}

/*
 * Source : MA in /u2/tm/src/srcgeo/julian_lib/ 
 */
/*
 * Inputs : year_number : year number (4 digits) julian_day : integer day number or
 * julian day (1..366) 
 */
/*
 * Outputs : day_of_month : day of the month (1..31) month_number : month number (1..12) 
 */
/*
 * The procedure "julian_to_date" does the reverse operation of the procedure
 * "make_julian_day" i.e. computes the month number and the respective day of month from 
 * the information on year and integer day number. Returns 0 if OK, 1 otherwise. 
 */
int sg1_day_of_year_to_ymd(int year, int day_of_year, int *month,
        int *day_of_month)
{
    if (year <= 0)
        return 1;
    if (day_of_year < 1)
        return 1;

    int const * day_of_year_offset;

    if (is_leap_year(year)) {
        day_of_year_offset = SG1_MONTHLY_DAY_OF_YEAR_OFFSET_LEAP_YEAR;
    } else {
        day_of_year_offset = SG1_MONTHLY_DAY_OF_YEAR_OFFSET;
    }

    if (day_of_year > day_of_year_offset[12])
        return 1;

    int m = 12;
    while(day_of_year <= day_of_year_offset[--m]);

    *month = m+1;
    *day_of_month = day_of_year - day_of_year_offset[m];
    return 0;

}

/*
 * Source : 
 */
/*
 * Inputs : year_number : year number (4 digits) month_number : month number (1..12) 
 */
/*
 * Outputs : number_days_month : number of days in a month 
 */
/*
 * The procedure "nbdays_month" gives the number of days in a month, useful for monthly
 * calculations. Returns 0 if OK, 1 otherwise. 
 */
int sg1_nbdays_month(int year, int month, int *number_days_of_month)
{
    if (year <= 0)
        return 1;
    if ((month < 1) || (month > 12))
        return 1;

    *number_days_of_month = SG1_DAYS_PER_MONTH[month - 1];
    if (is_leap_year(year))
        *number_days_of_month += 1;
    return 0;
}

/*
 * Source : 
 */
/*
 * Inputs : month_number : month number (1..12) month_name : name of month (3 characters
 * only, jan..dec) 
 */
/*
 * Outputs : month_name : name of the month abbreviated with 3 characters (jan..dec) 
 */
/*
 * The procedure "number_to_name_month" converts the month number into the corresponding 
 * month name. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_number_to_name_month (int month_number, char *month_name)
{
   char tab_name[][4] =
    { "   ", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
    "sep", "oct", "nov", "dec"
  };
  int ier;

  ier = 1;
  if ((month_number > 0) && (month_number < 13))
    {
      ier = 0;
      sprintf (month_name, "%s", tab_name[month_number]);
    }

  return (ier);
}

/**
 * The procedure "Day_Angle" expresses the integer day number as an angle (in radians)
 * from 12:00 hours on the day 31st December. A year length of 365.2422 days is used.
 * Returns 0 if OK, 1 otherwise.
 *
 * @input day_of_year: the day number within the year in [1,366]
 * @return Day_Angle
 **/
double sg1_day_angle(int day_of_year) {
    return day_of_year * 2.0 * SG1_PI_LOW_PRECISION / 365.2422;
}

/**
 * The procedure "Day_Angle" expresses the integer day number as an angle (in radians)
 * from 12:00 hours on the day 31st December. A year length of 365.2422 days is used.
 * Returns 0 if OK, 1 otherwise.
 *
 * @input day_of_year: the day number within the year in [1,366]
 * @return Day_Angle
 **/
int sg1_day_angle(int day_of_year, double *day_angle)
{
    if ((day_of_year < 1) || (day_of_year > 366))
        return 1;
    *day_angle = sg1_day_angle(day_of_year);
    return 0;
}

/**
 * The procedure "declination_sun" computes the solar declination at noon in solar time
 * (in radians). A single (average) value per day -at noon- is adequate for pratical
 * calculations. The noon declination depends on longitude, as noon occurs earlier if
 * longitude is East of Greenwi²ch, and later if it is West. The chosen algorithm uses
 * 1957 as base year; it is basically a truncated Fourier series with six harmonics.
 * Returns 0 if OK, 1 otherwise.
 *
 * Sources : Bourges, B., 1985. Improvement in solar declination computation. Solar
 * Energy, 35 (4), 367-369. Carvalho, M.J. and Bourges, B., 1986. Program Eufrad 2.0 -
 * User's Guide. Project EUFRAT final scientific report, Contract EN3S-0111-F, Solar
 * Energy and Development in the European Community, pp. 12.1-12.74. Duffie, J.A. and
 * Beckman, W.A., 1980. Solar Engineering of Thermal Processes. Wiley-Interscience, New
 * York.
 *
 * @input year: the year number usualy 4 digits
 * @input day_of_year: the number of the day within year
 * @input lambda: the lingitude of the cite in radians
 * @return declination of the sun in radians
 **/
double sg1_declination_sun(int year, int day_of_year, double lambda)
{
    double const b1 = 0.0064979;
    double const b2 = 0.4059059;
    double const b3 = 0.0020054;
    double const b4 = -0.0029880;
    double const b5 = -0.0132296;
    double const b6 = 0.0063809;
    double const b7 = 0.0003508;

    /*
     * n0 : spring-equinox time expressed in days from the beginning of the year i.e.
     * the time in decimal days elapsing from 00:00 hours Jan 1st to the spring equinox
     * at Greenwich in a given year
     * t1 : time in days, from the spring equinox
     * 0.5 represents the decimal day number at noon on Jan 1st at Greenwich
     */
    double n0 = 78.8946 + 0.2422 * (year - 1957) - floor(0.25 * (year - 1957));
    double t1 = -0.5 - lambda / (2 * SG1_PI_LOW_PRECISION) - n0;
    double w0 = 2 * SG1_PI_LOW_PRECISION / 365.2422;
    double wt = w0 * (day_of_year + t1);

    return b1 + b2 * sin(wt) + b3 * sin(2 * wt) + b4 * sin(3 * wt)
            + b5 * cos(wt) + b6 * cos(2 * wt) + b7 * cos(3 * wt);
}

/**
 * The procedure "declination_sun" computes the solar declination at noon in solar time
 * (in radians). A single (average) value per day -at noon- is adequate for pratical
 * calculations. The noon declination depends on longitude, as noon occurs earlier if
 * longitude is East of Greenwi²ch, and later if it is West. The chosen algorithm uses
 * 1957 as base year; it is basically a truncated Fourier series with six harmonics.
 * Returns 0 if OK, 1 otherwise.
 *
 * Sources : Bourges, B., 1985. Improvement in solar declination computation. Solar
 * Energy, 35 (4), 367-369. Carvalho, M.J. and Bourges, B., 1986. Program Eufrad 2.0 -
 * User's Guide. Project EUFRAT final scientific report, Contract EN3S-0111-F, Solar
 * Energy and Development in the European Community, pp. 12.1-12.74. Duffie, J.A. and
 * Beckman, W.A., 1980. Solar Engineering of Thermal Processes. Wiley-Interscience, New
 * York.
 *
 * @input year: the year number usualy 4 digits
 * @input day_of_year: the number of the day within year
 * @input lambda: the lingitude of the cite in radians
 * @return declination of the sun in radians
 **/
int sg1_declination_sun(int year, int day_of_year, double lambda,
        double *delta)
{
    if ((day_of_year < 1) || (day_of_year > 366))
        return 1;
  *delta = sg1_declination_sun(year, day_of_year, lambda);
  return 0;
}


/*
 * Source : Gruter (Ed.) (1984) 
 */
/*
 * Inputs : month_number : month number (1..12) 
 */
/*
 * Outputs : delta_month : solar declination angle (in radians) 
 */
/*
 * The procedure "declination_sun_month" computes the noon solar declination (in
 * radians) in solar time, with a simplified form, in two cases: type_use=0 : for
 * estimating monthly mean global solar radiation type_use=1 : for estimating monthly
 * mean maximum global solar radiation The integer day number to be selected in each case 
 * for the computations is given by two tables. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_declination_sun_month (int month_number, int type_use, double *delta_month)
{
  const double deg_rad = (SG1_PI_LOW_PRECISION / 180.0);	/* converts decimal degrees into radians */
  int tab_julian_day[12] =
    { 17, 46, 75, 105, 135, 162, 198, 228, 259, 289, 319, 345 };
  int tab_julian_day_max[12] =
    { 29, 57, 89, 119, 150, 173, 186, 217, 248, 278, 309, 339 };
  int ier, julian_day;
  double day_angle, jm, c1, c2, c3, c4;

  ier = 1;
  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 0)
	julian_day = tab_julian_day[month_number - 1];
      if (type_use == 1)
	julian_day = tab_julian_day_max[month_number - 1];
      ier = 0;
    }

  ier = sg1_day_angle (julian_day, &day_angle);
  if (ier != 0)
    return (ier);

  jm = day_angle;
  c1 = 0.3978;
  c2 = 80.2 * deg_rad;		/* 1.4000 in SSA manual */
  c3 = 1.92 * deg_rad;		/* 0.0355 in SSA manual */
  c4 = 2.80 * deg_rad;		/* 0.0489 in SSA manual */

  *delta_month = asin (c1 * sin (jm - c2 + c3 * sin (jm - c4)));
  return (ier);
}


/*
 * Source : 
 */
/*
 * Inputs : t : solar time i.e. LAT (0..24 decimal hours) 
 */
/*
 * Outputs : omega : solar hour angle (in radians) 
 */
/*
 * The procedure "solar_hour_angle" supplies the solar hour angle (in radians). By
 * convention the hour angle is negative before noon and positive after noon Returns 0 if 
 * OK, 1 otherwise. 
 */
 int
sg1_solar_hour_angle (double t, double *omega)
{
  int ier;

  ier = 1;
  if ((t >= 0.0) && (t <= 24.0))
    {
      ier = 0;
      *omega = (t - 12.0) * SG1_PI_LOW_PRECISION / 12.0;
    }
  return (ier);
}


/*
 * Source : 
 */
/*
 * Inputs : omega : solar hour angle (in radians) 
 */
/*
 * Outputs : t : solar time i.e. LAT (0..24 decimal hours) 
 */
/*
 * The procedure "omega_to_LAT" does the reverse operation of the procedure
 * "solar_hour_angle" i.e. computes the solar time (in decimal hours) from the solar
 * hour angle (in radians). Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_omega_to_LAT (double omega, double *t)
{
  int ier;

  ier = 1;

  if ((omega >= -SG1_PI_LOW_PRECISION) && (omega <= SG1_PI_LOW_PRECISION))
    {
      ier = 0;
      *t = 12.0 * (1.0 + omega / SG1_PI_LOW_PRECISION);
    }
  return (ier);
}

 double
sg1_geogr_to_geoce (double phi_g)
{
  double phi;
  double CC = 0.99330552;	/* Correction factor for converting geographic */
  /*
   * into geocentric latitude. CC=(Rpole/Requator)**2 
   */
  /*
   * Rpole=6356.752, Requator=6378.137 
   */
  if ((phi_g >= -(SG1_PI_LOW_PRECISION / 2.0 - 0.0002)) || (phi_g <= (SG1_PI_LOW_PRECISION / 2.0 - 0.0002)))
    phi = atan (tan (phi_g) * CC);
  else
    phi = phi_g;
  return (phi);
}

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude (in radians, positive to North) delta : solar declination
 * angle (in radians) t : solar time i.e. LAT (0..24 decimal hours) 
 */
/*
 * Outputs : omega : solar hour angle (in radians) 
 */
/*
 * The procedure "solar_hour_angle_h" supplies an average value of the solar hour angle
 * (in radians) for a whole solar hour, taking into account only the portion of the solar 
 * hour with the sun standing above the horizon. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_solar_hour_angle_h (double phi_g, double delta, double t, double *omega)
{
  int ier;
  double omega_sr, omega_ss, omega1, omega2;

  ier = 1;
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);

  omega1 = (t - 1.0 - 12.0) * SG1_PI_LOW_PRECISION / 12.0;
  if (omega1 < omega_sr)
    omega1 = omega_sr;

  omega2 = (t - 12.0) * SG1_PI_LOW_PRECISION / 12.0;
  if (omega2 > omega_ss)
    omega2 = omega_ss;

  *omega = (omega1 + omega2) / 2.0;

  return (ier);
}

/*********************************/
/*
 * SUNRISE, SUNSET AND DAYLENGTH 
 */
/*********************************/


/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) delta : solar
 * declination angle (in radians) gamma_riset : solar elevation near sunrise/sunset: -
 * set to 0.0 for astronomical sunrise/sunset - set to -1.0 for refraction corrected
 * sunrise/sunset. 
 */
/*
 * Outputs : omega_sr : sunrise solar hour angle (in radians) omega_ss : sunset solar
 * hour angle (in radians) 
 */
/*
 * The procedure "sunrise_hour_angle" supplies the sunrise and sunset hour angles (in
 * radians). Due to the dimension of the solar disk and the effect of the atmospheric
 * refraction, the edge of the solar disk will just appear (disappear) at the horizon at 
 * sunrise (at sunset) when the calculated astronomical elevation is 50'. Returns 0 if
 * OK, 1 otherwise. 
 */
 int
sg1_sunrise_hour_angle (double phi_g, double delta, double gamma_riset,
		    double *omega_sr, double *omega_ss)
{
   double deg_rad = (SG1_PI_LOW_PRECISION / 180.0);	/* converts decimal degrees into radians */
  int ier;
  double horizon, max_delta, cos_omegas, omegas;
  double phi;

  ier = 1;
  if ((gamma_riset == 0.0) || (gamma_riset == -1.0))
    ier = 0;
  horizon = (-50.0 / 60.0) * deg_rad;	/* horizon, -50' in radians */
  if (gamma_riset >= horizon)
    horizon = gamma_riset;

  phi = sg1_geogr_to_geoce (phi_g);
  max_delta = 23.45 * deg_rad;

  if ((fabs (phi) < (SG1_PI_LOW_PRECISION / 2.0)) && (fabs (delta) <= max_delta) && (ier == 0))
    {
      cos_omegas =
	(sin (horizon) - (sin (phi) * sin (delta))) / (cos (phi) * cos (delta));
      ier = 0;
    }
  else
    ier = 1;

  if (fabs (cos_omegas) < 1.0)
    omegas = acos (cos_omegas);
  if (cos_omegas >= 1.0)	/* the sun is always below the horizon : polar night */
    omegas = 0.0;
  if (cos_omegas <= -1.0)	/* the sun is always above the horizon : polar day */
    omegas = SG1_PI_LOW_PRECISION;

  *omega_sr = -omegas;
  *omega_ss = omegas;
  return (ier);
}


/**
 * @input phi: geocentric phi in radians (latitude)
 * @input delta: sun declination in radians
 * @return omega at sunset
 **/
double sg1_sunset(double phi, double delta)
{
  return atan2(-1.0 * sin (phi) * sin (delta), cos (phi) * cos (delta));
}


/*
 * Source : 
 */
/*
 * Inputs : omega_sr : sunrise hour angle (in radians) omega_ss : sunset hour angle (in
 * radians) 
 */
/*
 * Outputs : t_sr : time of astronomical sunrise (in decimal hours) t_ss : time of
 * astronomical sunset (in decimal hours) S0 : astronomical daylength (in decimal hours) 
 */
/*
 * The procedure "timerise_daylength" supplies the times of astronomical sunrise and
 * sunset, and the astronomical daylength, all in LAT decimal hours. Returns 0 if OK, 1
 * otherwise. 
 */
 int
sg1_timerise_daylength (double omega_sr, double omega_ss, double *t_sr,
		    double *t_ss, double *S0)
{
  int ier;

  ier = 1;
  if ((omega_sr >= -SG1_PI_LOW_PRECISION) && (omega_sr <= 0.0) && (omega_ss >= 0.0)
      && (omega_ss <= SG1_PI_LOW_PRECISION))
    {
      ier = 0;
      /*
       * alternative way
       */
      /*
       * ier = omega_to_LAT(omega_sr,&t_sr); if(ier == 0) ier ==
       * omega_to_LAT(omega_ss,&t_ss); if(ier != 0) return(ier);
       */
      *t_sr = 12.0 + omega_sr * 12.0 / SG1_PI_LOW_PRECISION;
      *t_ss = 12.0 + omega_ss * 12.0 / SG1_PI_LOW_PRECISION;
      *S0 = *t_ss - *t_sr;
    }
  return (ier);
}

/****************************/
/*
 * CHANGING THE TIME SYSTEM 
 */
/****************************/

/*
 * Source : Gruter (ed.) (1984) 
 */
/*
 * Inputs : day_angle : day angle (in radians) lambda : longitude of the site (in
 * radians, positive to East) lambda_ref : reference longitude of the time zone (in
 * radians) summer_corr : correction for summer time (integer hours) 
 */
/*
 * Outputs : dt : Offset between local mean time (LMT) and local apparent time (LAT) (in
 * decimal hours) 
 */
/*
 * The procedure "LMT_to_LAT computes the difference (in decimal hours) between the LAT
 * (local apparent time) and the LMT (local mean time or clock time) systems at solar
 * noon. Two stages: - the first stage calculates the equation of time, ET, wich allows
 * for perturbations in the rotational and angular orbital speed of the Earth. - the
 * second stage handles the difference between the longitude of the site under
 * consideration and the reference time zone longitude for the site. A summer time
 * correction must be added for some countries. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_LMT_to_LAT (double day_angle, double lambda, double lambda_ref, int summer_corr,
	    double *dt)
{
  const double deg_rad = (SG1_PI_LOW_PRECISION / 180.0);	/* converts decimal degrees into radians */
  int ier;
  double a1, a2, a3, a4, ET;

  ier = 1;
  a1 = -0.128;
  a2 = -0.165;
  a3 = 2.80 * deg_rad;
  a4 = 19.70 * deg_rad;
  if ((day_angle > 0.0) && (day_angle < (2.0 * SG1_PI_LOW_PRECISION * 1.0021)) &&
      (fabs (lambda) <= SG1_PI_LOW_PRECISION) && (fabs (lambda_ref) <= SG1_PI_LOW_PRECISION))
    {
      ier = 0;
      ET = a1 * sin (day_angle - a3) + a2 * sin (2.0 * day_angle + a4);
      *dt = ET + ((lambda - lambda_ref) * 12.0 / SG1_PI_LOW_PRECISION) - (double) summer_corr;
    }

  return (ier);
}

/*
 * Source : 
 */
/*
 * Inputs : UT : Universal Time (in decimal hours) day_angle : day angle (in radians)
 * lambda : longitude of the site (in radians, positive to East) 
 */
/*
 * Outputs : LAT : local apparent time or solar time or true solar time (TST) (in decimal
 * hours) 
 */
/*
 * The procedure "UT_to_LAT computes the conversion of the UT (Universal time) into the
 * LAT (local apparent time) systems at solar noon (in decimal hours). First, the
 * equation of time, ET, is computed (in decimal hours), wich allows for perturbations in 
 * the rotational and angular orbital speed of the Earth. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_UT_to_LAT (double UT, double day_angle, double lambda, double *LAT)
{
  const double deg_rad = (SG1_PI_LOW_PRECISION / 180.0);	/* converts decimal degrees into radians */
  int ier;
  double a1, a2, a3, a4, ET;

  ier = 1;
  a1 = -0.128;
  a2 = -0.165;
  a3 = 2.80 * deg_rad;
  a4 = 19.70 * deg_rad;
  if ((UT >= 0.0) && (UT <= 24.0) && (day_angle > 0.0) &&
      (day_angle < (2.0 * SG1_PI_LOW_PRECISION * 1.0021)) && (fabs (lambda) <= SG1_PI_LOW_PRECISION))
    {
      ier = 0;
      ET = a1 * sin (day_angle - a3) + a2 * sin (2.0 * day_angle + a4);
      *LAT = UT + ET + (lambda * 12.0 / SG1_PI_LOW_PRECISION);
      if (*LAT < 0)
	*LAT += 24.0;
      if (*LAT > 24.0)
	*LAT -= 24.0;
    }

  return (ier);
}

/**********************************/
/*
 * POSITION OF THE SUN IN THE SKY 
 */
/**********************************/

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) delta : solar
 * declination angle (in radians) omega : solar hour angle (in radians) 
 */
/*
 * Outputs : gamma : solar altitude angle (in radians) theta : solar zenithal angle (in
 * radians) 
 */
/*
 * The procedure "elevation_zenith_sun" computes the solar elevation (or altitude) angle 
 * and the solar zenithal (or incidence) angle. These two angles are complementary.
 * Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_elevation_zenith_sun (double phi_g, double delta, double omega, double *gamma,
		      double *theta)
{
  int ier;
  double omega_sr, omega_ss;
  double phi;

  ier = 1;
  phi = sg1_geogr_to_geoce (phi_g);
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  if ((omega < omega_sr) || (omega > omega_ss))
    *gamma = 0.0;
  else
    *gamma =
      asin (sin (phi) * sin (delta) + cos (phi) * cos (delta) * cos (omega));
  if (*gamma < 0.0)
    *gamma = 0.0;

  *theta = (SG1_PI_LOW_PRECISION / 2.0) - *gamma;

  return (ier);
}

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) delta : solar
 * declination angle (in radians) omega : solar hour angle (in radians) gamma : solar
 * altitude angle (in radians) 
 */
/*
 * Outputs : alpha : solar azimuthal angle (in radians) 
 */
/*
 * The procedure "azimuth_sun" computes the solar azimuth angle in the Northern
 * hemisphere. The azimuth angle has a positive value when the sun is to the west of
 * South, i.e. during the afternoon in solar time. For the Southern hemisphere, the
 * azimuth angle is measured from North. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_azimuth_sun (double phi_g, double delta, double omega, double gamma,
	     double *alpha)
{
  int ier;
  double cos_as, sin_as, x;
  double phi;

  ier = 0;
  phi = sg1_geogr_to_geoce (phi_g);
  cos_as = (sin (phi) * sin (gamma) - sin (delta)) / (cos (phi) * cos (gamma));
  if (phi < 0.0)
    cos_as = -cos_as;		/* Southern hemisphere */
  sin_as = cos (delta) * sin (omega) / cos (gamma);
  
  if (cos_as > 1.0) 
	  cos_as = 1.0;
  
  if (cos_as < -1.0) 
	  cos_as = -1.0;

  x = acos (cos_as);
  if (fabs (x) > SG1_PI_LOW_PRECISION)
    ier = 1;

  if (sin_as >= 0.0)
    *alpha = x;
  else
    *alpha = -x;

  return (ier);
}

/********************************/
/*
 * EXTRATERRESTRIAL IRRADIATION 
 */
/********************************/

/*
 * Source : Gruter (ed.) (1984) 
 */
/*
 * Inputs : day_angle : day angle (in radians) 
 */
/*
 * Outputs : eccentricity : correction for Earth orbit eccentricity 
 */
/*
 * The procedure "corr_distance" computes the correction for the variation of sun-earth
 * distance from its mean value (also known as eccentricity). It is a fucntion of time,
 * but a single (average) value per day is enough for practical calculations. Returns 0
 * if OK, 1 otherwise. 
 */
 int
sg1_corr_distance (double day_angle, double *eccentricity)
{
  const double deg_rad = (SG1_PI_LOW_PRECISION / 180.0);	/* converts decimal degrees into radians */
  int ier;
  double a;

  ier = 1;
  a = 2.80 * deg_rad;
  if ((day_angle >= 0.0) && (day_angle <= (2.0 * SG1_PI_LOW_PRECISION * 1.0021)))
    {
      ier = 0;
      *eccentricity = 1.0 + 0.03344 * cos (day_angle - a);
    }
  return (ier);
}


/*
 * Source : 
 */
/*
 * Inputs : IOj : extraterrestrial solar irradiance normal to beam for day j; I0j=I0*fj
 * theta : solar incidence angle or solar zenithal angle 
 */
/*
 * Outputs : G0 : extraterrestrial global solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "G0_normal" delivers the extraterrestrial solar irradiance normal to
 * beam for day j. Returns 0 if OK, 1 otherwise. 
 */
 int
sg1_G0_normal (double I0j, double theta, double *G0)
{
  *G0 = I0j * cos (theta);
  return 0;
}

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) eccentricity :
 * correction for Earth orbit eccentricity delta : solar declination angle (in radians)
 * omega1 : solar hour angle at beginning of the period (in radians) omega2 : solar hour
 * angle at end of the period (in radians) 
 */
/*
 * Outputs : G0_12 : extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "G0_general" delivers the extraterrestrial solar irradiation incident
 * on an horizontal surface in the general case (in Wh/m2). Returns 0 if OK, 1 otherwise 
 */
 int
sg1_G0_general (double phi_g, double eccentricity, double delta,
	    double omega1, double omega2, double *G0_12)
{
  int ier;
  double omega_sr, omega_ss, a, b1, b2, c;
  double phi;

  ier = 1;
  phi = sg1_geogr_to_geoce (phi_g);
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  if (omega1 < omega_sr)
    omega1 = omega_sr;
  if (omega2 < omega_sr)
    omega2 = omega_sr;
  if (omega1 > omega_ss)
    omega1 = omega_ss;
  if (omega2 > omega_ss)
    omega2 = omega_ss;

  if (omega2 <= omega1)
    *G0_12 = 0.0;
  else
    {
      a = SG1_I0 * eccentricity * SG1_DAY_LENGTH / (2.0 * SG1_PI_LOW_PRECISION);
      b1 = sin (phi) * sin (delta) * (omega2 - omega1);
      b2 = cos (phi) * cos (delta) * (sin (omega2) - sin (omega1));
      c = a * (b1 + b2);
      if (c < 0.0)
	*G0_12 = 0.0;
      else
	*G0_12 = c;
    }

  return (ier);
}

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) eccentricity :
 * correction for Earth orbit eccentricity delta : solar declination angle (in radians) 
 */
/*
 * Outputs : G0d : daily extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "G0_day" delivers the extraterrestrial solar irradiation incident on an 
 * horizontal surface in case of daily values (in Wh/m2), i.e. omega1 = omega_sr =
 * -omega_ss et omega2 = omega_ss. Returns 0 if OK, 1 otherwise. REMARK: It is a special 
 * case of G0_general with the sunrise and sunset angles as integration limits. 
 */
 int
sg1_G0_day (double phi_g, double eccentricity, double delta, double *G0d)
{
  int ier;
  double omega_sr, omega_ss, a, b;
  double phi;

  ier = 1;
  phi = sg1_geogr_to_geoce (phi_g);
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = SG1_I0 * eccentricity * SG1_DAY_LENGTH / SG1_PI_LOW_PRECISION;
  /*
   * b = cos(phi) * cos(delta) * (sin(omega_ss) - omega_ss * cos(omega_ss)); 
   */
  b =
    sin (phi) * sin (delta) * omega_ss +
    cos (phi) * cos (delta) * sin (omega_ss);
  *G0d = a * b;

  return (ier);
}

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) eccentricity :
 * correction for Earth orbit eccentricity delta : solar declination (in radians) 
 */
/*
 * Outputs : G0h[1..24] : 24 hourly extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "G0_hours_profile" delivers the extraterrestrial solar irradiation
 * incident on an horizontal surface in case of hourly values, for the 24 integral hours
 * in a given day (in Wh/m2), i.e. |omega1-omega2| = Pi/12. Returns 0 if OK, 1 otherwise 
 */
 int
sg1_G0_hours_profile (double phi_g, double eccentricity, double delta, double *G0h)
{
  int ier, i;
  double omega_sr, omega_ss, a, b1, b2;
  double phi;
  double t1, t2, omega1, omega2;

  ier = 1;
  phi = sg1_geogr_to_geoce (phi_g);
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = SG1_I0 * eccentricity * SG1_DAY_LENGTH / (2.0 * SG1_PI_LOW_PRECISION);
  b1 = sin (phi) * sin (delta);
  b2 = cos (phi) * cos (delta);

  for (i = 0; i < 24; i++)
    {
      t1 = (double) (i + 1) - 1.0;
      ier = sg1_solar_hour_angle (t1, &omega1);
      if (ier != 0)
	return (ier);
      t2 = (double) (i + 1);
      ier = sg1_solar_hour_angle (t2, &omega2);
      if (ier != 0)
	return (ier);

      if ((omega2 < omega_sr) || (omega1 > omega_ss))
	G0h[i] = 0.0;
      else
	{
	  if (omega1 < omega_sr)
	    omega1 = omega_sr;
	  if (omega2 > omega_ss)
	    omega2 = omega_ss;
	  G0h[i] =
	    a * (b1 * (omega2 - omega1) + b2 * (sin (omega2) - sin (omega1)));
	}
    }

  return (ier);
}

/*
 * Source : 
 */
/*
 * Inputs : phi_g : latitude of site (in radians, positive to North) eccentricity :
 * correction for Earth orbit eccentricity delta : solar declination (in radians) t :
 * solar time i.e. LAT (0..24 decimal hours) 
 */
/*
 * Outputs : G0h : hourly extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "G0_hour" delivers the extraterrestrial solar irradiation incident on
 * an horizontal surface for a specific hour in a given day (in Wh/m2), i.e.
 * |omega1-omega2| = Pi/12. t is taken as the mid hour for computation of the hourly
 * value of extraterrestrial solar irradiation. Returns 0 if OK, 1 otherwise 
 */
 int
sg1_G0_hour (double phi_g, double eccentricity, double delta, double t, double *G0h)
{
  int ier;
  double omega_sr, omega_ss, a, b1, b2;
  double t1, t2, omega1, omega2;
  double phi;

  ier = 1;
  phi = sg1_geogr_to_geoce (phi_g);
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = SG1_I0 * eccentricity * SG1_DAY_LENGTH / (2.0 * SG1_PI_LOW_PRECISION);
  b1 = sin (phi) * sin (delta);
  b2 = cos (phi) * cos (delta);

  t1 = t - 1.0;
  ier = sg1_solar_hour_angle (t1, &omega1);
  if (ier != 0)
    return (ier);
  t2 = t;
  ier = sg1_solar_hour_angle (t2, &omega2);
  if (ier != 0)
    return (ier);

  if (omega1 < omega_sr)
    omega1 = omega_sr;
  if (omega2 < omega_sr)
    omega2 = omega_sr;

  if (omega2 <= omega1)
    *G0h = 0.0;
  else
    {
      *G0h = a * (b1 * (omega2 - omega1) + b2 * (sin (omega2) - sin (omega1)));
      if (*G0h < 0.0)
	*G0h = 0.0;
    }

  return (ier);
}

/***********************************************/
/*
 * MONTHLY AVERAGES OF SOLAR INPUTS PARAMETERS 
 */
/***********************************************/

/*
 * Source : 
 */
/*
 * Inputs : month_number : month number (1..12) year_number : year number (4 digits)
 * phi_g : latitude of site (in radians, positive to North) lambda : longitude of site
 * (in radians, positive to East) gamma_riset : solar elevation near sunrise/sunset: -
 * set to 0.0 for astronomical sunrise/sunset - set to -1.0 for refraction corrected
 * sunrise/sunset. 
 */
/*
 * Outputs : monthly average of... day_angle_m : ... day angle (in radians) delta_m : ... 
 * solar declination angle (in radians) omega_ss_m : ... sunset hour angle (in radians)
 * S0_m : ... astronomical daylength (in decimal hours) eccentricity_m : ... eccentricity
 * G0d_m : ... daily extraterrestrial irradiation (in Wh/m2) G0h_m[1..24] : ... 24 hourly
 * extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "monthly_averages" computes directly the monthly average values of
 * solar parameters : day angle (in radians), eccentricity, declination (in radians),
 * sunset hour angle (in radians), daylength (in decimal hours), daily extraterrestrial
 * irradiation (in Wh/m2) and the 24 hourly extraterrestrial solar irradiation (in
 * Wh/m2) . Returns 0 if OK, 1 otherwise 
 */
 int
sg1_monthly_averages (int month_number, int year_number,
		  double phi_g, double lambda, double gamma_riset,
		  double *day_angle_m, double *delta_m, double *omega_ss_m,
		  double *S0_m, double *eccentricity_m, double *G0d_m,
		  double *G0h_m)
{
  int ier, i, day_of_month, number_days_month, julian_day;
  double day_angle, delta, omega_sr, omega_ss, t_sr, t_ss, S0, eccentricity,
    G0d;
  double G0h[24];
  double nbd_m;

  /*
   * Initialization 
   */
  *day_angle_m = 0.0;
  *delta_m = 0.0;
  *omega_ss_m = 0.0;
  *S0_m = 0.0;
  *eccentricity_m = 0.0;
  *G0d_m = 0.0;
  for (i = 0; i < 24; i++)
    G0h_m[i] = 0.0;

  ier = 1;
  ier = sg1_nbdays_month (year_number, month_number, &number_days_month);
  if (ier != 0)
    return (ier);

  for (day_of_month = 1; day_of_month <= number_days_month; day_of_month++)
    {
      ier =
	sg1_ymd_to_day_of_year (day_of_month, month_number, year_number, &julian_day);
      if (ier == 0)
	ier = sg1_day_angle (julian_day, &day_angle);
      if (ier == 0)
	ier = sg1_declination_sun (year_number, julian_day, lambda, &delta);
      if (ier == 0)
	ier =
	  sg1_sunrise_hour_angle (phi_g, delta, gamma_riset, &omega_sr, &omega_ss);
      if (ier == 0)
	ier = sg1_timerise_daylength (omega_sr, omega_ss, &t_sr, &t_ss, &S0);
      if (ier == 0)
	ier = sg1_corr_distance (day_angle, &eccentricity);
      if (ier == 0)
	ier = sg1_G0_day (phi_g, eccentricity, delta, &G0d);
      if (ier == 0)
	ier = sg1_G0_hours_profile (phi_g, eccentricity, delta, G0h);
      if (ier != 0)
	return (ier);

      /*
       * OR 
       */
      /*
       * ier =
       * solar_parameters_day(day_of_month,month_number,year_number,phi_g,lambda,gamma_riset,&day_angle,&delta,&omega_ss,&S0,&eccentricity,&G0d,G0h); 
       */
      if (ier != 0)
	return (ier);

      /*
       * REMARK: In the original procedure test on G0d: if(*G0d > 0) 
       */
      *day_angle_m = *day_angle_m + day_angle;
      *delta_m = *delta_m + delta;
      *omega_ss_m = *omega_ss_m + omega_ss;
      *S0_m = *S0_m + S0;
      *eccentricity_m = *eccentricity_m + eccentricity;
      *G0d_m = *G0d_m + G0d;
      for (i = 0; i < 24; i++)
	G0h_m[i] = G0h_m[i] + G0h[i];
    }

  nbd_m = (double) number_days_month;
  *day_angle_m = *day_angle_m / nbd_m;
  *delta_m = *delta_m / nbd_m;
  *omega_ss_m = *omega_ss_m / nbd_m;
  *S0_m = *S0_m / nbd_m;
  *eccentricity_m = *eccentricity_m / nbd_m;
  *G0d_m = *G0d_m / nbd_m;
  for (i = 0; i < 24; i++)
    G0h_m[i] = G0h_m[i] / nbd_m;

  return (ier);
}

/****************************************************/
/*
 * YEARLY AVERAGES OF MONTHLY SOLAR PARAMETERS 
 */
/*
 * (LONG TERM MEANS OF MONTHLY MEANS OF DAILY SUMS) 
 */
/****************************************************/

/*
 * Source : 
 */
/*
 * Inputs : month_number : month number (1..12) year_start : starting year of the
 * considered period (4 digits) year_end : ending year of the considered period (4
 * digits) phi_g : latitude of site (in radians, positive to North) lambda : longitude of 
 * site (in radians, positive to East) gamma_riset : solar elevation near sunrise/sunset:
 * - set to 0.0 for astronomical sunrise/sunset - set to -1.0 for refraction corrected
 * sunrise/sunset. 
 */
/*
 * Outputs : yearly average of... day_angle_y : ... day angle (in radians) delta_y : ...
 * solar declination angle (in radians) omega_ss_y : ... sunset hour angle (in radians)
 * S0_y : ... astronomical daylength (in decimal hours) eccentricity_y : ... eccentricity
 * G0d_y : ... daily extraterrestrial irradiation (in Wh/m2) G0h_y[1..24] : ... 24 hourly
 * extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "yearly_averages" computes directly the yearly average over a defined
 * period of years of monthly average values of solar parameters : day angle (in
 * radians), eccentricity, declination (in radians), sunset hour angle (in radians),
 * daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2) and 24
 * hourly extraterrestrial solar irradiation (in Wh/m2). Returns 0 if OK, 1 otherwise 
 */
 int
sg1_yearly_averages (int month_number, int year_start, int year_end,
		 double phi_g, double lambda, double gamma_riset,
		 double *day_angle_y, double *delta_y, double *omega_ss_y,
		 double *S0_y, double *eccentricity_y, double *G0d_y,
		 double *G0h_y)
{
  int ier, i, year_number;
  double day_angle_m, delta_m, omega_ss_m, S0_m, eccentricity_m, G0d_m,
    G0h_m[25];
  double number_of_years;

  number_of_years = (double) (year_end - year_start + 1.);
  ier = 1;

  /*
   * Initialization of parameters 
   */
  *day_angle_y = 0.0;
  *delta_y = 0.0;
  *omega_ss_y = 0.0;
  *S0_y = 0.0;
  *eccentricity_y = 0.0;
  *G0d_y = 0.0;
  for (i = 0; i < 24; i++)
    G0h_y[i] = 0.0;

  for (year_number = year_start; year_number <= year_end; year_number++)
    {
      ier =
	sg1_monthly_averages (month_number, year_number, phi_g, lambda,
			  gamma_riset, &day_angle_m, &delta_m,
			  &omega_ss_m, &S0_m, &eccentricity_m, &G0d_m, G0h_m);
      if (ier != 0)
	return (ier);

      *day_angle_y = *day_angle_y + day_angle_m;
      *delta_y = *delta_y + delta_m;
      *omega_ss_y = *omega_ss_y + omega_ss_m;
      *S0_y = *S0_y + S0_m;
      *eccentricity_y = *eccentricity_y + eccentricity_m;
      *G0d_y = *G0d_y + G0d_m;
      for (i = 0; i < 24; i++)
	G0h_y[i] = G0h_y[i] + G0h_m[i];
      /*
       * printf("year_number= %4d G0d_m = %8.2f (Wh/m2) G0d_y = %8.2f
       * (Wh/m2)\n",year_number,G0d_m,*G0d_y); for(i=0;i<24;i++) printf("hh= %2d
       * (hours)\tG0h_m = %8.2f (Wh/m2)\tG0h_y = %8.2f
       * (Wh/m2)\n",i+1,G0h_m[i],G0h_y[i]); 
       */
    }

  *day_angle_y = *day_angle_y / number_of_years;
  *delta_y = *delta_y / number_of_years;
  *omega_ss_y = *omega_ss_y / number_of_years;
  *S0_y = *S0_y / number_of_years;
  *eccentricity_y = *eccentricity_y / number_of_years;
  *G0d_y = *G0d_y / number_of_years;
  for (i = 0; i < 24; i++)
    G0h_y[i] = G0h_y[i] / number_of_years;

  return (ier);
}

/*********************************************/
/*
 * SOLAR INPUTS PARAMETERS FOR A CERTAIN DAY 
 */
/*********************************************/

/*
 * Source : 
 */
/*
 * Inputs : day_of_month : day of the month (1..31) month_number : month number (1..12)
 * year_number : year number (4 digits) phi_g : latitude of site (in radians, positive to 
 * North) lambda : longitude of site (in radians, positive to East) gamma_riset : solar
 * elevation near sunrise/sunset: - set to 0.0 for astronomical sunrise/sunset - set to
 * -1.0 for refraction corrected sunrise/sunset. 
 */
/*
 * Outputs : day_angle : day angle (in radians) delta : solar declination angle (in
 * radians) omega_ss : sunset hour angle (in radians) S0 : astronomical daylength (in
 * decimal hours) eccentricity : eccentricity G0d : daily extraterrestrial irradiation
 * (in Wh/m2) G0h[1..24] : 24 hourly extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "solar_parameters_day" computes the solar geometry related values for a 
 * certain day : day angle (in radians), eccentricity, declination (in radians), sunset
 * hour angle (in radians), daylength (in decimal hours), daily extraterrestrial
 * irradiation (in Wh/m2) and the 24 hourly extraterrestrial solar irradiation (in
 * Wh/m2). Returns 0 if OK, 1 otherwise. REMARK: gamma_riset set to 0.0 in the original 
 * procedure by Aguiar.
 */
 int
sg1_solar_parameters_day (int day_of_month, int month_number, int year_number,
		      double phi_g, double lambda, double gamma_riset,
		      double *day_angle, double *delta, double *omega_ss,
		      double *S0, double *eccentricity, double *G0d,
		      double *G0h)
{
  int ier, julian_day;
  double omega_sr, t_sr, t_ss;

  ier = 1;
  ier = sg1_ymd_to_day_of_year (day_of_month, month_number, year_number, &julian_day);
  if (ier == 0)
    ier = sg1_day_angle (julian_day, day_angle);
  if (ier == 0)
    ier = sg1_declination_sun (year_number, julian_day, lambda, delta);
  if (ier == 0)
    ier = sg1_sunrise_hour_angle (phi_g, *delta, gamma_riset, &omega_sr, omega_ss);
  if (ier == 0)
    ier = sg1_timerise_daylength (omega_sr, *omega_ss, &t_sr, &t_ss, S0);
  if (ier == 0)
    ier = sg1_corr_distance (*day_angle, eccentricity);
  if (ier == 0)
    ier = sg1_G0_day (phi_g, *eccentricity, *delta, G0d);
  if (ier == 0 && *G0d > 0.0)
    ier = sg1_G0_hours_profile (phi_g, *eccentricity, *delta, G0h);

  return (ier);
}

/******************************************************************/
/*
 * SOLAR INPUTS PARAMETERS FOR MONTHLY AVERAGE IRRADIATION MODELS 
 */
/******************************************************************/

/*
 * Source : 
 */
/*
 * Inputs : month_number : month number (1..12) phi_g : latitude of site (in radians,
 * positive to North) gamma_riset : solar elevation near sunrise/sunset: - set to 0.0 for 
 * astronomical sunrise/sunset - set to -1.0 for refraction corrected sunrise/sunset. 
 */
/*
 * Outputs : average ... for the given month day_angle_avg : day angle (in radians)
 * delta_avg : solar declination angle (in radians) omega_ss_avg : sunset hour angle (in
 * radians) S0_avg : astronomical daylength (in decimal hours) eccentricity_avg :
 * eccentricity G0d_avg : daily extraterrestrial irradiation (in Wh/m2) G0h_avg[1..24] :
 * 24 hourly extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "solar_parameters_acg" computes the solar geometry related values for
 * monthly average irradiation models : day angle (in radians), eccentricity,
 * declination (in radians), sunset hour angle (in radians), daylength (in decimal
 * hours), daily extraterrestrial irradiation (in Wh/m2) and the 24 hourly
 * extraterrestrial solar irradiation (in Wh/m2). Returns 0 if OK, 1 otherwise. REMARK: 
 * gamma_riset set to 0.0 in the original procedure by Aguiar.
 */
 int
sg1_solar_parameters_avg (int month_number,
		      double phi_g, double gamma_riset,
		      double *day_angle_avg, double *delta_avg,
		      double *omega_ss_avg, double *S0_avg,
		      double *eccentricity_avg, double *G0d_avg,
		      double *G0h_avg)
{
  int ier, julian_day;
  double omega_sr, t_sr, t_ss;

  /*
   * recommended values of day number for estimating monthly mean global solar
   * radiation 
   */
  int tab_julian_day[12] =
    { 17, 46, 75, 105, 135, 162, 198, 228, 259, 289, 319, 345 };
  int type_use;

  ier = 1;
  type_use = 0;			/* for estimating monthly mean global solar radiation */

  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 0)
	julian_day = tab_julian_day[month_number - 1];
      ier = 0;
    }
  if (ier == 0)
    ier = sg1_day_angle (julian_day, day_angle_avg);
  if (ier == 0)
    ier = sg1_declination_sun_month (month_number, type_use, delta_avg);
  if (ier == 0)
    ier =
      sg1_sunrise_hour_angle (phi_g, *delta_avg, gamma_riset, &omega_sr,
			  omega_ss_avg);
  if (ier == 0)
    ier = sg1_timerise_daylength (omega_sr, *omega_ss_avg, &t_sr, &t_ss, S0_avg);
  if (ier == 0)
    ier = sg1_corr_distance (*day_angle_avg, eccentricity_avg);
  if (ier == 0)
    ier = sg1_G0_day (phi_g, *eccentricity_avg, *delta_avg, G0d_avg);
  if (ier == 0)
    ier = sg1_G0_hours_profile (phi_g, *eccentricity_avg, *delta_avg, G0h_avg);

  return (ier);
}

/******************************************************************/
/*
 * SOLAR INPUTS PARAMETERS FOR MONTHLY MAXIMUM IRRADIATION MODELS 
 */
/******************************************************************/

/*
 * Source : 
 */
/*
 * Inputs : month_number : month number (1..12) phi_g : latitude of site (in radians,
 * positive to North) gamma_riset : solar elevation near sunrise/sunset: - set to 0.0 for 
 * astronomical sunrise/sunset - set to -1.0 for refraction corrected sunrise/sunset. 
 */
/*
 * Outputs : average ... for the given month day_angle_max : day angle (in radians)
 * delta_max : solar declination angle (in radians) omega_ss_max : sunset hour angle (in
 * radians) S0_max : astronomical daylength (in decimal hours) eccentricity_max :
 * eccentricity G0d_max : daily extraterrestrial irradiation (in Wh/m2) G0h_max[1..24] :
 * 24 hourly extraterrestrial solar irradiation (in Wh/m2) 
 */
/*
 * The procedure "solar_parameters_acg" computes the solar geometry related values for
 * monthly average irradiation models : day angle (in radians), eccentricity,
 * declination (in radians), sunset hour angle (in radians), daylength (in decimal
 * hours), daily extraterrestrial irradiation (in Wh/m2) and the 24 hourly
 * extraterrestrial solar irradiation (in Wh/m2). Returns 0 if OK, 1 otherwise. REMARK: 
 * gamma_riset set to 0.0 in the original procedure by Aguiar.
 */
 int
sg1_solar_parameters_max (int month_number,
		      double phi_g, double gamma_riset,
		      double *day_angle_max, double *delta_max,
		      double *omega_ss_max, double *S0_max,
		      double *eccentricity_max, double *G0d_max,
		      double *G0h_max)
{
  int ier, julian_day;
  double omega_sr, t_sr, t_ss;

  /*
   * recommended values of day number for estimating monthly mean maximum global solar 
   * radiation 
   */
  int tab_julian_day_max[12] =
    { 29, 57, 89, 119, 150, 173, 186, 217, 248, 278, 309, 339 };
  int type_use;

  ier = 1;
  type_use = 1;			/* for estimating monthly mean global solar radiation */

  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 1)
	julian_day = tab_julian_day_max[month_number - 1];
      ier = 0;
    }
  if (ier == 0)
    ier = sg1_day_angle (julian_day, day_angle_max);
  if (ier == 0)
    ier = sg1_declination_sun_month (month_number, type_use, delta_max);
  if (ier == 0)
    ier =
      sg1_sunrise_hour_angle (phi_g, *delta_max, gamma_riset, &omega_sr,
			  omega_ss_max);
  if (ier == 0)
    ier = sg1_timerise_daylength (omega_sr, *omega_ss_max, &t_sr, &t_ss, S0_max);
  if (ier == 0)
    ier = sg1_corr_distance (*day_angle_max, eccentricity_max);
  if (ier == 0)
    ier = sg1_G0_day (phi_g, *eccentricity_max, *delta_max, G0d_max);
  if (ier == 0)
    ier = sg1_G0_hours_profile (phi_g, *eccentricity_max, *delta_max, G0h_max);

  return (ier);
}


/*
 * Inputs : phi_g : geographic latitude of site (in radians, positive to North) delta :
 * delta : solar declination angle (radian) omega_ss : sunset hour angle (in radians)
 * beta : tilt angle of the inclined flat plane
 * alpha : azimuth of the inclined flat
 * plane 
 */
/*
 * Outputs : day_angle_max : day angle (in radians) delta_max : solar declination angle
 * (in radians) omega_ss_max : sunset hour angle (in radians) S0_max : astronomical
 * daylength (in decimal hours) eccentricity_max : eccentricity G0d_max : daily
 * extraterrestrial irradiation (in Wh/m2) G0h_max[1..24] : 24 hourly extraterrestrial
 * solar irradiation (in Wh/m2) 
 */
/*
 * Compute the solar angles omega when the Sun is visible (located in front of) by the
 * tilted surface [omega is related to time TST by: tTST = 12(1+omega/pi)] an area of 4
 * values is returned: [$omega1, $omega2, $omega3, $omega4] If all values are -999, then 
 * the sun is not visible by the plane at any time of the day If $omega3 and $omega4 are
 * -999, and $omega1 and $omega2 not, sun is visible during the period [$omega1;$omega2]
 * If all values are different from -999,sun is visible during the
 * periods:[$omega1;$omega2] and [$omega3;$omega4] 
 */
 int
sg1_intervals_omega_tilted_plane (double phi_g, double delta, double omega_ss,
			      double beta, double alpha, double *v_om,
			      int *p_nb)
{
  int ierr = 0;
  double A, B, C;
  double precision = 1e-4;
  double epsilon_z = 1e-6;
  double mAB, sAB;
  int nzA, nzB, nzC, cas;
  double d;
  double cosTheta_wm, cosTheta_wt;
  double wm, wt, wa, wb, wmin, wmax;
  double omega_sr;
  double r1, r2;
  double sin_wa, cos_wa;
  double sin_wb, cos_wb;
  double sin_w1, cos_w1;
  double sin_w2, cos_w2;
  int pa, pb;
  double phi;

  phi = sg1_geogr_to_geoce (phi_g);

  v_om[0] = -999.0;
  v_om[1] = -999.0;
  v_om[2] = -999.0;
  v_om[3] = -999.0;
  *p_nb = 0;

  if (fabs (omega_ss) < precision)
    {
      return ierr;
    }

  omega_sr = -omega_ss;

  A = cos (delta) * (cos (phi) * cos (beta) +
		     sin (phi) * sin (beta) * cos (alpha));
  
  B = cos (delta) * sin (beta) * sin (alpha);
  C = sin (delta) * (sin (phi) * cos (beta) -
		      cos (phi) * sin (beta) * cos (alpha));

  nzA = (int) (fabs (A) > epsilon_z);
  nzB = (int) (fabs (B) > epsilon_z);
  nzC = (int) (fabs (C) > epsilon_z);

  cas = nzC + nzB * 2 + nzA * 4;


  switch (cas)
    {
    case 0:
      /*
       * ANGLES = [omega_sr omega_ss]; 
       */
      v_om[0] = omega_sr;
      v_om[1] = omega_ss;
      break;
    case 1:
      if (C > 0)
	{
	  /*
	   * ANGLES = [omega_sr omega_ss]; 
	   */
	  v_om[0] = omega_sr;
	  v_om[1] = omega_ss;
	}
      break;
    case 2:
      if (B > 0)
	{
	  /*
	   * ANGLES = [0 omega_ss]; 
	   */
	  v_om[0] = 0.0;
	  v_om[1] = omega_ss;
	}
      else
	{
	  /*
	   * ANGLES = [omega_sr 0]; 
	   */
	  v_om[0] = omega_sr;
	  v_om[1] = 0.0;
	}
      break;
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:

      mAB = A * A + B * B;
      sAB = sqrt (mAB);

      d = fabs (C) / sAB;
      if (fabs (d - 1) < precision)
	{

	  if (C > 0.0)
	    {
	      wa = -asin (B);
	    }
	  else
	    {
	      wa = asin (B);
	    }

	  if ((wa < omega_sr + precision) | (wa > omega_ss - precision))
	    {
	      wt = (omega_sr + omega_ss) / 2.0;
	      cosTheta_wt = cos (wt) * A + sin (wt) * B + C;
	      if (cosTheta_wt > precision)
		{
		  /*
		   * ANGLES = [omega_sr omega_ss]; 
		   */
		  v_om[0] = omega_sr;
		  v_om[1] = omega_ss;
		}
	      else if (fabs (cosTheta_wt) < precision)
		{
		  ierr = 1;
		}
	    }
	  else
	    {
	      /*
	       * Solution soit [omega_sr wa] ou [wa omega_ss] �d�inir Test au
	       * centre et au bord omega_sr ou omega_ss du plus grand intervalle 
	       */
	      if (omega_ss - wa > wa - omega_sr)
		{
		  wt = (omega_ss + wa) / 2;
		  wm = omega_ss;
		}
	      else
		{
		  wt = (omega_sr + wa) / 2;
		  wm = omega_sr;
		}
	      cosTheta_wm = cos (wm) * A + sin (wm) * B + C;
	      cosTheta_wt = cos (wt) * A + sin (wt) * B + C;

	      if (fabs (cosTheta_wm) > fabs (cosTheta_wt))
		{
		  cosTheta_wt = cosTheta_wm;
		}
	      if (cosTheta_wt > precision)
		{
		  /*
		   * ANGLES = [wa omega_ss]; 
		   */
		  v_om[0] = wa;
		  v_om[1] = omega_ss;
		}
	      else if (cosTheta_wt < -precision)
		{
		  /*
		   * ANGLES = [omega_sr wa]; 
		   */
		  v_om[0] = omega_sr;
		  v_om[1] = wa;
		}
	      else
		{
		  ierr = 1;
		}
	    }
	}
      else if (d >= 1 + precision)
	{

	  wt = (omega_sr + omega_ss) / 2.0;
	  cosTheta_wt = cos (wt) * A + sin (wt) * B + C;
	  if (cosTheta_wt > precision)
	    {
	      /*
	       * ANGLES = [omega_sr omega_ss]; 
	       */
	      v_om[0] = omega_sr;
	      v_om[1] = omega_ss;
	    }
	  else if (fabs (cosTheta_wt) < precision)
	    {
	      ierr = 1;
	    }
	}
      else
	{
	  if (cas != 6)
	    {

	      r1 = fabs (A) * sqrt (mAB - C * C);
	      sin_w1 = -(B * C + r1) / mAB;
	      sin_w2 = -(B * C - r1) / mAB;

	      r2 = fabs (B) * sqrt (mAB - C * C);
	      cos_w1 = -(A * C + r2) / mAB;
	      cos_w2 = -(A * C - r2) / mAB;

	      if (fabs (sin_w1 * sin_w1 + cos_w1 * cos_w1 - 1.0) > epsilon_z)
		{
		  sin_wa = sin_w1;
		  cos_wa = cos_w2;
		  sin_wb = sin_w2;
		  cos_wb = cos_w1;
		}
	      else
		{
		  sin_wa = sin_w1;
		  cos_wa = cos_w1;
		  sin_wb = sin_w2;
		  cos_wb = cos_w2;
		}

	      wa = atan2 (sin_wa, cos_wa);
	      wb = atan2 (sin_wb, cos_wb);
	    }
	  else
	    {
	      wa = atan (-A / B);
	      wb = SG1_PI_LOW_PRECISION + wa;
	      if (wb > SG1_PI_LOW_PRECISION)
		{
		  wb = wb - 2.0 * SG1_PI_LOW_PRECISION;
		}
	    }

	  wmin = MIN (wa, wb);
	  wmax = MAX (wa, wb);
	  wa = wmin;
	  wb = wmax;

	  pa = 0;
	  if (wa < omega_sr + precision)
	    {
	      pa = -1;
	    }
	  else if (wa > omega_ss - precision)
	    {
	      pa = 1;
	    }

	  pb = 0;
	  if (wb < omega_sr + precision)
	    {
	      pb = -1;
	    }
	  else if (wb > omega_ss - precision)
	    {
	      pb = 1;
	    }

	  if ((pa != 0) && (pb != 0))
	    {
	      /*
	       * [wa omega_sr omega_ss wb] 
	       */
	      wt = (omega_sr + omega_ss) / 2.0;
	      cosTheta_wt = cos (wt) * A + sin (wt) * B + C;
	      if (cosTheta_wt > precision)
		{
		  /*
		   * ANGLES = [omega_sr omega_ss]; 
		   */
		  v_om[0] = omega_sr;
		  v_om[1] = omega_ss;
		}
	      else if (fabs (cosTheta_wt) < precision)
		{
		  ierr = 1;
		}
	    }
	  else if ((pa == 0) && (pb == 0))
	    {
	      /*
	       * [omega_sr wa wb omega_ss] 
	       */
	      wt = (wa + wb) / 2.0;
	      cosTheta_wt = cos (wt) * A + sin (wt) * B + C;
	      if (cosTheta_wt > precision)
		{
		  /*
		   * ANGLES = [wa wb]; 
		   */
		  v_om[0] = wa;
		  v_om[1] = wb;
		}
	      else if (cosTheta_wt < -precision)
		{
		  /*
		   * ANGLES = [omega_sr wa wb omega_ss]; 
		   */
		  v_om[0] = omega_sr;
		  v_om[1] = wa;
		  v_om[2] = wb;
		  v_om[3] = omega_ss;
		}
	      else
		{
		  ierr = 1;
		}
	    }
	  else if (pa == 0)
	    {
	      /*
	       * [omega_sr wa omega_ss wb] 
	       */
	      wt = (omega_sr + wa) / 2.0;
	      cosTheta_wt = cos (wt) * A + sin (wt) * B + C;
	      if (cosTheta_wt > precision)
		{
		  /*
		   * ANGLES = [omega_sr wa]; 
		   */
		  v_om[0] = omega_sr;
		  v_om[1] = wa;
		}
	      else if (cosTheta_wt < -precision)
		{
		  /*
		   * ANGLES = [wa omega_ss]; 
		   */
		  v_om[0] = wa;
		  v_om[1] = omega_ss;
		}
	      else
		{
		  ierr = 1;
		}
	    }
	  else
	    {
	      /*
	       * [wa omega_sr wb omega_ss] 
	       */
	      wt = (omega_sr + wb) / 2.0;
	      cosTheta_wt = cos (wt) * A + sin (wt) * B + C;
	      if (cosTheta_wt > precision)
		{
		  /*
		   * ANGLES = [omega_sr wb]; 
		   */
		  v_om[0] = omega_sr;
		  v_om[1] = wb;
		}
	      else if (cosTheta_wt < -precision)
		{
		  /*
		   * ANGLES = [wb omega_ss]; 
		   */
		  v_om[0] = wb;
		  v_om[1] = omega_ss;
		}
	      else
		{
		  ierr = 1;
		}
	    }
	}
    }

  if (v_om[0] == -999.0)
    {
      *p_nb = 0;
    }
  else if (v_om[2] == -999.0)
    {
      *p_nb = 2;
    }
  else
    {
      *p_nb = 0;
    }

  /* Patch pour faire fonctionner le code de Mireille */
  v_om[0] = omega_sr;
  v_om[1] = omega_ss;
  *p_nb = 2;

  return (ierr);

}


double sg1_ymd_to_julian_day(int year, int month, int day_of_month)
{
       if (month <= 2) {
               month = month + 12;
               year = year - 1;
       }
       return 1721028.0 + day_of_month + floor((153 * month - 2) / 5) + 365 * year + floor(year / 4)
                   - floor(year / 100) + floor(year / 400) + 12.0 / 24.0 - 0.5;
}



/***************************************/
/* POSITION OF THE SUN IN THE SKY FAST */
/***************************************/

void sg1_init_solar_geometry_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double phi_g, double delta)
{
	double phi;

	p_sgf->phi_g = phi_g;
	p_sgf->delta = delta;

	phi = sg1_geogr_to_geoce(phi_g);
	p_sgf->phi = phi;

	p_sgf->sin_phi = sin(phi);
	p_sgf->cos_phi = cos(phi);
	p_sgf->sin_delta = sin(delta);
	p_sgf->cos_delta = cos(delta);

	p_sgf->sin_phi_sin_delta = p_sgf->sin_phi*p_sgf->sin_delta;
	p_sgf->cos_phi_cos_delta = p_sgf->cos_phi*p_sgf->cos_delta;

}

void sg1_deftilt_solar_geometry_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double alpha, double beta)
{
	p_sgf->alpha = alpha;
	p_sgf->sin_alpha = sin(alpha);
	p_sgf->cos_alpha = cos(alpha);
	p_sgf->beta = beta;
	p_sgf->sin_beta = sin(beta);
	p_sgf->cos_beta = cos(beta);

	if (p_sgf->phi >= 0.0) {
		p_sgf->A = p_sgf->cos_delta * (p_sgf->cos_phi * p_sgf->cos_beta + p_sgf->sin_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
		p_sgf->B = p_sgf->cos_delta * p_sgf->sin_beta * p_sgf->sin_alpha;
		p_sgf->C = p_sgf->sin_delta * (p_sgf->sin_phi * p_sgf->cos_beta - p_sgf->cos_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
	} else {
		p_sgf->A = p_sgf->cos_delta * (p_sgf->cos_phi * p_sgf->cos_beta - p_sgf->sin_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
		p_sgf->B = p_sgf->cos_delta * p_sgf->sin_beta * p_sgf->sin_alpha;
		p_sgf->C = p_sgf->sin_delta * (p_sgf->sin_phi * p_sgf->cos_beta + p_sgf->cos_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
	}
}

void sg1_elevation_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega, double *p_gamma)
{
	*p_gamma = asin(p_sgf->sin_phi_sin_delta + p_sgf->cos_phi_cos_delta*cos_omega);
}

void sg1_elevation_zenith_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega, double *p_gamma, double *p_theta)
{
	*p_gamma = asin(p_sgf->sin_phi_sin_delta + p_sgf->cos_phi_cos_delta*cos_omega);
	*p_theta = SG1_PI_2 - *p_gamma;
}

void sg1_azimuth_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double sin_omega, double gamma, double *p_alpha)
{
	double cos_as, sin_as, x;
	double cos_gamma = cos(gamma);

	cos_as = (p_sgf->sin_phi * sin(gamma) - p_sgf->sin_delta) / (p_sgf->cos_phi * cos_gamma);
	if (p_sgf->phi < 0.0)
		cos_as = -cos_as; /* Southern hemisphere */

	sin_as = p_sgf->cos_delta * sin_omega / cos_gamma;

	if (cos_as > 1.0)
		cos_as = 1.0;
	else if (cos_as < -1.0)
		cos_as = -1.0;

	x = acos(cos_as);
	if (sin_as >= 0.0)
		*p_alpha = x;
	else
		*p_alpha = -x;

}

void sg1_cos_incident_angle_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega, double sin_omega, double *p_costhetai)
{
	*p_costhetai = p_sgf->A*cos_omega + p_sgf->B*sin_omega + p_sgf->C;
}

