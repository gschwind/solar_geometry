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

#include <cmath>
#include <algorithm>



namespace sg1 {

/* convert degree to radians, using SG1_PI_LOW_PRECISION */
inline static constexpr double RAD_LOW_PRECISSION(double a)
{
	return a * sg1::PI_LOW_PRECISION / 180.0;
}

/* convert radians to degree, using SG1_PI_LOW_PRECISION */
inline static constexpr double DEG_LOW_PRECISSION(double a)
{
	return a * 180.0 / sg1::PI_LOW_PRECISION;
}

/* the last value is used for sg1_day_of_year_to_ymd */
static int const MONTHLY_DAY_OF_YEAR_OFFSET[13] = {
        0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
};

/* the last value is used for sg1_day_of_year_to_ymd */
static int const MONTHLY_DAY_OF_YEAR_OFFSET_LEAP_YEAR[13] = {
        0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366
};

static int DAYS_PER_MONTH[12] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

/*
 * recommended values of day number for estimating monthly mean global solar
 * radiation
 */
static int const DAY_OF_YEAR_FOR_MONTHLY_MEAN_ESTIMATION[12] = {
		17, 46, 75, 105, 135, 162, 198, 228, 259, 289, 319, 345
};

/*
 * recommended values of day number for estimating monthly mean maximum global solar
 * radiation
 */
static int const DAY_OF_YEAR_FOR_MONTHLY_MAX_ESTIMATION[12] = {
		29, 57, 89, 119, 150, 173, 186, 217, 248, 278, 309, 339
};


static constexpr double const DELTA_MAX = RAD_LOW_PRECISSION(23.45);


inline static int is_leap_year(int const year) {
    return (((year % 4) == 0) && ((year % 100) != 0)) || ((year % 400) == 0);
}

int nbdays_month(int year, int month)
{
    int number_days_of_month = DAYS_PER_MONTH[month - 1];
    if (is_leap_year(year)&&(month==2))
        number_days_of_month += 1;
    return number_days_of_month;
}


double day_angle(int day_of_year) {
    return day_of_year * 2.0 * sg1::PI_LOW_PRECISION / 365.2422;
}


double declination_sun(int year, int day_of_year, double lambda)
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
    double t1 = -0.5 - lambda / (2 * sg1::PI_LOW_PRECISION) - n0;
    double w0 = 2 * sg1::PI_LOW_PRECISION / 365.2422;
    double wt = w0 * (day_of_year + t1);

    return b1 + b2 * sin(wt) + b3 * sin(2 * wt) + b4 * sin(3 * wt)
            + b5 * cos(wt) + b6 * cos(2 * wt) + b7 * cos(3 * wt);
}


double solar_hour_angle(double t)
{
    return (t - 12.0) * sg1::PI_LOW_PRECISION / 12.0;
}


/**
 * Supplies the solar time in hours in [0,24].
 *
 * @input omega: solar_hour_angle in radians
 * @return solar time i.e. LAT (0..24 decimal hours)
 **/
double omega_to_LAT(double omega)
{
    return 12.0 + omega * 12.0 / sg1::PI_LOW_PRECISION;
}


double sunset(double phi, double delta)
{
  return acos(- tan(phi) * tan(delta));
}


std::tuple<int,int,int> julian_day_to_ymd(int jd)
{
	double H, L, N, I, J, K;

	L = jd + 68569.0;
	N = floor(4 * L / 146097.0);
	L = L - floor((146097.0 * N + 3.0) / 4.0);
	I = floor(4000 * (L + 1) / 1461001.0);
	L = L - floor(1461.0 * I / 4.0) + 31.0;

	J = floor(80.0 * L / 2447.0);
	K = L - floor(2447.0 * J / 80.0);
	L = floor(J / 11.0);
	J = J + 2.0 - 12.0 * L;
	I = 100.0 * (N - 49.0) + I + L;

	return std::tuple<int,int,int>(I,J,K);

}


double gamma_sun(double phi, double delta, double omega)
{
    return asin (sin (phi) * sin (delta) + cos (phi) * cos (delta) * cos (omega));
}


double corr_distance(double day_angle)
{
    double const a = RAD_LOW_PRECISSION(2.80);
    return 1.0 + 0.03344 * cos(day_angle - a);
}


inline static double compute_ET(double day_angle)
{
    double const a1 = -0.128;
    double const a2 = -0.165;
    double const a3 = sg1::RAD_LOW_PRECISSION(2.80);
    double const a4 = sg1::RAD_LOW_PRECISSION(19.70);
    return a1 * sin (day_angle - a3) + a2 * sin (2.0 * day_angle + a4);
}


double azimuth_sun(double phi, double delta, double omega, double gamma)
{
    double cos_as = sin(phi) * sin(gamma) - sin(delta);
    double sin_as = cos(delta) * sin(omega) * cos(phi);
    /* Check hemisphere and swap the azimuth if needed */
    if (phi < 0.0)
        cos_as = -cos_as;
    return atan2(sin_as, cos_as);
}


double geogr_to_geoce(double phi_g)
{
    /**
     * Correction factor for converting geographic
     * into geocentric latitude. CC=(Rpole/Requator)**2
     * with Rpole=6356.752, Requator=6378.137
     **/
    double const CC = 0.99330552;
    return atan2(sin(phi_g)*CC, cos(phi_g));
}

} // namespace sg1

int sg1_ymd_to_day_of_year(int year, int month, int day_of_month,
        int * const day_of_year)
{
    int julien;

    if (year <= 0)
        return 1;
    if ((month < 1) || (month > 12))
        return 1;

    /* Technically this is not required */
    if ((day_of_month < 1) || (day_of_month > 31))
        return 1;

    *day_of_year = day_of_month + sg1::MONTHLY_DAY_OF_YEAR_OFFSET[month - 1];
    if (sg1::is_leap_year(year) && (month > 2))
        *day_of_year += 1;
    return 0;
}


int sg1_day_of_year_to_ymd(int year, int day_of_year, int *month,
        int *day_of_month)
{
    if (year <= 0)
        return 1;
    if (day_of_year < 1)
        return 1;

    int const * day_of_year_offset;

    if (sg1::is_leap_year(year)) {
        day_of_year_offset = sg1::MONTHLY_DAY_OF_YEAR_OFFSET_LEAP_YEAR;
    } else {
        day_of_year_offset = sg1::MONTHLY_DAY_OF_YEAR_OFFSET;
    }

    if (day_of_year > day_of_year_offset[12])
        return 1;

    int m = 12;
    while(day_of_year <= day_of_year_offset[--m]);

    *month = m+1;
    *day_of_month = day_of_year - day_of_year_offset[m];
    return 0;

}

int sg1_nbdays_month(int year, int month, int *number_days_of_month)
{
    if (year <= 0)
        return 1;
    if ((month < 1) || (month > 12))
        return 1;

    *number_days_of_month = sg1::nbdays_month(year, month);
    return 0;
}


int sg1_day_angle(int day_of_year, double *day_angle)
{
    if ((day_of_year < 1) || (day_of_year > 366))
        return 1;
    *day_angle = sg1::day_angle(day_of_year);
    return 0;
}


int sg1_declination_sun(int year, int day_of_year, double lambda,
        double *delta)
{
    if ((day_of_year < 1) || (day_of_year > 366))
        return 1;
  *delta = sg1::declination_sun(year, day_of_year, lambda);
  return 0;
}

static inline double _declination_sun_month(double day_angle)
{
	double const c1 = 0.3978;
	double const c2 = sg1::RAD_LOW_PRECISSION(80.2); /* 1.4000 in SSA manual */
	double const c3 = sg1::RAD_LOW_PRECISSION(1.92); /* 0.0355 in SSA manual */
	double const c4 = sg1::RAD_LOW_PRECISSION(2.80); /* 0.0489 in SSA manual */

	return asin(c1 * sin(day_angle - c2 + c3 * sin(day_angle - c4)));
}

int sg1_declination_sun_month_avg(int month_number, double *delta_month)
{
	double day_of_year =
			sg1::DAY_OF_YEAR_FOR_MONTHLY_MEAN_ESTIMATION[month_number - 1];

	double day_angle;

	if (sg1_day_angle(day_of_year, &day_angle) != 0)
		return 1;

	*delta_month = _declination_sun_month(day_angle);
	return 0;
}


int sg1_declination_sun_month_max(int month_number, double *delta_month)
{
	double day_of_year =
			sg1::DAY_OF_YEAR_FOR_MONTHLY_MAX_ESTIMATION[month_number - 1];

	double day_angle;

	if (sg1_day_angle(day_of_year, &day_angle) != 0)
		return 1;

	*delta_month = _declination_sun_month(day_angle);
	return 0;
}

int sg1_solar_hour_angle(double t, double *omega)
{
    if ((t < 0.0) || (t > 24.0))
        return 1;
    *omega = sg1::solar_hour_angle(t);
    return 0;
}


int sg1_omega_to_LAT(double omega, double *t)
{
    if (fabs(omega) > sg1::PI_LOW_PRECISION)
        return 1;
    *t = sg1::omega_to_LAT(omega);
    return 0;
}


int sg1_solar_hour_angle_h(double phi_g, double delta, double t,
        double *omega)
{
  int ier;
  double omega_sr, omega_ss, omega1, omega2;

  ier = 1;
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);

  omega1 = sg1::solar_hour_angle(t - 1.0);
  if (omega1 < omega_sr)
    omega1 = omega_sr;

  omega2 = sg1::solar_hour_angle(t);
  if (omega2 > omega_ss)
    omega2 = omega_ss;

  *omega = (omega1 + omega2) / 2.0;

  return (ier);
}


int sg1_sunrise_hour_angle(double phi_g, double delta, double gamma_riset,
        double *omega_sr, double *omega_ss)
{
  int ier;
  double horizon, cos_omegas, omegas;
  double phi;

  ier = 1;
  if ((gamma_riset == 0.0) || (gamma_riset == -1.0))
    ier = 0;
  horizon = sg1::RAD_LOW_PRECISSION(-50.0 / 60.0);	/* horizon, -50' in radians */
  if (gamma_riset >= horizon)
    horizon = gamma_riset;

  phi = sg1_geogr_to_geoce (phi_g);

  if ((fabs (phi) < (sg1::PI_LOW_PRECISION / 2.0)) && (fabs (delta) <= sg1::DELTA_MAX) && (ier == 0))
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
    omegas = sg1::PI_LOW_PRECISION;

  *omega_sr = -omegas;
  *omega_ss = omegas;
  return (ier);
}



int sg1_timerise_daylength(double omega_sr, double omega_ss, double *t_sr,
        double *t_ss, double *S0)
{
    if (omega_sr < -sg1::PI_LOW_PRECISION || omega_sr > 0.0)
        return 1;
    if (omega_ss > sg1::PI_LOW_PRECISION || omega_ss < 0.0)
        return 1;
    *t_sr = sg1::omega_to_LAT(omega_sr);
    *t_ss = sg1::omega_to_LAT(omega_ss);
    *S0 = *t_ss - *t_sr;
    return 0;
}


int sg1_LMT_to_LAT(double day_angle, double lambda, double lambda_ref,
        int summer_corr, double *dt)
{
    if (day_angle <= 0.0 || day_angle >= (2.0 * sg1::PI_LOW_PRECISION * 1.0021))
        return 1;
    if (fabs(lambda) > sg1::PI_LOW_PRECISION)
        return 1;
    if (fabs(lambda_ref) > sg1::PI_LOW_PRECISION)
        return 1;
    double ET = sg1::compute_ET(day_angle);
    *dt = ET + ((lambda - lambda_ref) * 12.0 / sg1::PI_LOW_PRECISION)
            - summer_corr;
    return 0;
}


int sg1_UT_to_LAT(double UT, double day_angle, double lambda, double *LAT)
{
    if (day_angle <= 0.0 || day_angle >= (2.0 * sg1::PI_LOW_PRECISION * 1.0021))
        return 1;
    if ((UT < 0.0) || (UT > 24.0))
        return 1;
    if (fabs(lambda) > sg1::PI_LOW_PRECISION)
        return 1;

    double ET = sg1::compute_ET(day_angle);
    *LAT = UT + ET + (lambda * 12.0 / sg1::PI_LOW_PRECISION);
    if (*LAT < 0.0)
        *LAT += 24.0;
    if (*LAT > 24.0)
        *LAT -= 24.0;
    return 0;
}




int sg1_elevation_zenith_sun(double phi_g, double delta, double omega,
        double *gamma, double *theta)
{
    if (fabs(delta) > sg1::DELTA_MAX)
        return 1;

    double phi = sg1_geogr_to_geoce(phi_g);
    if (fabs(phi) >= (sg1::PI_LOW_PRECISION / 2.0))
        return 1;

    double omega_ss = sg1_sunset(phi, delta);

    if (fabs(omega) >= omega_ss) {
        *gamma = 0.0;
    } else {
        *gamma = sg1::gamma_sun(phi, delta, omega);
        *gamma = std::max(0.0, *gamma);
    }

    *theta = (sg1::PI_LOW_PRECISION / 2.0) - *gamma;

    return 0;
}


int sg1_azimuth_sun(double phi_g, double delta, double omega, double gamma,
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
  if (fabs (x) > sg1::PI_LOW_PRECISION)
    ier = 1;

  if (sin_as >= 0.0)
    *alpha = x;
  else
    *alpha = -x;

  return (ier);
}

int sg1_corr_distance(double day_angle, double *eccentricity)
{
    if ((day_angle < 0.0) || (day_angle > (2.0 * sg1::PI_LOW_PRECISION * 1.0021)))
        return 1;
    *eccentricity = sg1::corr_distance(day_angle);
    return 0;
}


int sg1_G0_normal(double I0j, double theta, double *G0)
{
  *G0 = I0j * cos (theta);
  return 0;
}


int sg1_G0_general(double phi_g, double eccentricity, double delta,
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
      a = sg1::I0 * eccentricity * sg1::DAY_LENGTH / (2.0 * sg1::PI_LOW_PRECISION);
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


int sg1_G0_day(double phi_g, double eccentricity, double delta, double *G0d)
{
  int ier;
  double omega_sr, omega_ss, a, b;
  double phi;

  ier = 1;
  phi = sg1_geogr_to_geoce (phi_g);
  ier = sg1_sunrise_hour_angle (phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = sg1::I0 * eccentricity * sg1::DAY_LENGTH / sg1::PI_LOW_PRECISION;
  /*
   * b = cos(phi) * cos(delta) * (sin(omega_ss) - omega_ss * cos(omega_ss)); 
   */
  b =
    sin (phi) * sin (delta) * omega_ss +
    cos (phi) * cos (delta) * sin (omega_ss);
  *G0d = a * b;

  return (ier);
}


int sg1_G0_hours_profile(double phi_g, double eccentricity, double delta,
        double *G0h)
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
  a = sg1::I0 * eccentricity * sg1::DAY_LENGTH / (2.0 * sg1::PI_LOW_PRECISION);
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


int sg1_G0_hour(double phi_g, double eccentricity, double delta, double t,
        double *G0h)
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
  a = sg1::I0 * eccentricity * sg1::DAY_LENGTH / (2.0 * sg1::PI_LOW_PRECISION);
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


int sg1_monthly_averages(int month_number, int year_number, double phi_g,
        double lambda, double gamma_riset, double *day_angle_m, double *delta_m,
        double *omega_ss_m, double *S0_m, double *eccentricity_m, double *G0d_m,
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


int sg1_yearly_averages(int month_number, int year_start, int year_end,
        double phi_g, double lambda, double gamma_riset, double *day_angle_y,
        double *delta_y, double *omega_ss_y, double *S0_y,
        double *eccentricity_y, double *G0d_y, double *G0h_y)
{
  int ier, i, year_number;
  double day_angle_m, delta_m, omega_ss_m, S0_m, eccentricity_m, G0d_m,
    G0h_m[24];
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


int sg1_solar_parameters_day(int day_of_month, int month_number,
        int year_number, double phi_g, double lambda, double gamma_riset,
        double *day_angle, double *delta, double *omega_ss, double *S0,
        double *eccentricity, double *G0d, double *G0h)
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


int sg1_solar_parameters_avg(int month_number, double phi_g, double gamma_riset,
		double *day_angle_avg, double *delta_avg, double *omega_ss_avg,
		double *S0_avg, double *eccentricity_avg, double *G0d_avg,
		double *G0h_avg)
{
	double omega_sr, t_sr, t_ss;


	int day_of_year =
			sg1::DAY_OF_YEAR_FOR_MONTHLY_MEAN_ESTIMATION[month_number - 1];

	int ier = 0;
	if (ier == 0)
		ier = sg1_day_angle(day_of_year, day_angle_avg);
	if (ier == 0)
		ier = sg1_declination_sun_month_avg(month_number, delta_avg);
	if (ier == 0)
		ier = sg1_sunrise_hour_angle(phi_g, *delta_avg, gamma_riset, &omega_sr,
				omega_ss_avg);
	if (ier == 0)
		ier = sg1_timerise_daylength(omega_sr, *omega_ss_avg, &t_sr, &t_ss,
				S0_avg);
	if (ier == 0)
		ier = sg1_corr_distance(*day_angle_avg, eccentricity_avg);
	if (ier == 0)
		ier = sg1_G0_day(phi_g, *eccentricity_avg, *delta_avg, G0d_avg);
	if (ier == 0)
		ier = sg1_G0_hours_profile(phi_g, *eccentricity_avg, *delta_avg,
				G0h_avg);

	return ier;
}


int sg1_solar_parameters_max(int month_number, double phi_g, double gamma_riset,
        double *day_angle_max, double *delta_max, double *omega_ss_max,
        double *S0_max, double *eccentricity_max, double *G0d_max,
        double *G0h_max)
{
  int ier, julian_day;
  double omega_sr, t_sr, t_ss;


  int type_use;

  ier = 1;
  type_use = 1;			/* for estimating monthly mean global solar radiation */

  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 1)
	julian_day = sg1::DAY_OF_YEAR_FOR_MONTHLY_MAX_ESTIMATION[month_number - 1];
      ier = 0;
    }
  if (ier == 0)
    ier = sg1_day_angle (julian_day, day_angle_max);
  if (ier == 0)
    ier = sg1_declination_sun_month_max(month_number, delta_max);
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


int sg1_intervals_omega_tilted_plane(double phi_g, double delta,
        double omega_ss, double beta, double alpha, double *v_om, int *p_nb)
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
	      wb = sg1::PI_LOW_PRECISION + wa;
	      if (wb > sg1::PI_LOW_PRECISION)
		{
		  wb = wb - 2.0 * sg1::PI_LOW_PRECISION;
		}
	    }

	  wmin = std::min(wa, wb);
	  wmax = std::max(wa, wb);
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

/* return the julian day at 12h */
int sg1_ymd_to_julian_day(int year, int month, int day_of_month)
{
	int k;
	double Y, M, D, H;

	Y = year;
	M = month;
	D = day_of_month;
	if (M < 3) {
		M += 12;
		Y -= 1;
	}

	return 1721028.0 + D + floor((153.0 * M - 2.0) / 5.0) + 365.0 * Y
			+ floor(Y / 4.0) - floor(Y / 100.0) + floor(Y / 400.0);

}

void sg1_julian_day_to_ymd(int jd, int * year, int * month, int * day_of_month)
{
	std::tie(*year,*month,*day_of_month) = sg1::julian_day_to_ymd(jd);
}



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
		p_sgf->A = p_sgf->cos_delta * (p_sgf->cos_phi * p_sgf->cos_beta
		         + p_sgf->sin_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
		p_sgf->B = p_sgf->cos_delta * p_sgf->sin_beta * p_sgf->sin_alpha;
		p_sgf->C = p_sgf->sin_delta * (p_sgf->sin_phi * p_sgf->cos_beta
		         - p_sgf->cos_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
	} else {
		p_sgf->A = p_sgf->cos_delta * (p_sgf->cos_phi * p_sgf->cos_beta
		         - p_sgf->sin_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
		p_sgf->B = p_sgf->cos_delta * p_sgf->sin_beta * p_sgf->sin_alpha;
		p_sgf->C = p_sgf->sin_delta * (p_sgf->sin_phi * p_sgf->cos_beta
		         + p_sgf->cos_phi * p_sgf->sin_beta * p_sgf->cos_alpha);
	}
}

void sg1_elevation_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega, double *p_gamma)
{
	*p_gamma = asin(p_sgf->sin_phi_sin_delta + p_sgf->cos_phi_cos_delta*cos_omega);
}

void sg1_elevation_zenith_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega, double *p_gamma, double *p_theta)
{
	*p_gamma = asin(p_sgf->sin_phi_sin_delta + p_sgf->cos_phi_cos_delta*cos_omega);
	*p_theta = sg1::PI_2 - *p_gamma;
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

