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

#ifndef __H_solar_geometry
#define __H_solar_geometry

#ifdef	__cplusplus
extern "C" {
#endif

/**
 * NOTATIONS
 * =========
 *
 * phi_g: geographic latitude of the site, positive to North
 * phi: geocentric latitude of the site, positive to North
 * lambda: longitude of the site, positive to East
 * delta: solar declination angle
 * omega: solar hour angle
 * gamma: solar altitude (or elevation) angle
 * theta: solar incidence (or zenithal) angle (ESRA --> zeta)
 * alpha: solar azimuthal angle (or psi)
 * t: solar time = true solar time (TST) = local apparent time (TST)
 * TST: local apparent time or solar time or true solar time (TST)
 *       --> this system of time offers the advantage of symmetry of the solar
 *       geometry about the north-south line
 * LMT: local mean time or clock time
 * UT: Universal Time, is GMT measured from Greenwich mean midnight
 * omega_sr: sunrise hour angle
 * omega_ss: sunset hour angle
 * t_sr: time of astronomical sunrise
 * t_ss: time of astronomical sunset
 * omega1: solar hour angle at beginning of the time period
 * omega2: solar hour angle at end of the time period
 * S0: astronomical daylength or astronomical sunshine duration
 * I0: solar constant = annual mean value of extraterrestrial direct solar
 *     irradiance G0 (1367.0 W/m2)
 * G0: extraterrestrial global solar irradiation (on an horizontal plane) =B0
 * G0h: hourly extraterrestrial solar irradiation (on an horizontal plane)
 * G0d: daily extraterrestrial solar irradiation (on an horizontal plane)
 *
 * Remarks:
 *  - All angles are computed in radians as a standard.
 *  - The basic trigonometric calculations on the position of the sun are
 *    carried out in TST.
 **/

/*********************************************
 *                                           *
 * G E O M E T R Y  O F  S O L A R  B E A M, *
 *                                           *
 * A  N E A R  P O I N T  S O U R C E        *
 *                                           *
 *********************************************/

/********************
 * BASIC PARAMETERS *
 ********************/


/**
 * The procedure "make_julian_day" converts a day given in day, month and year
 * into a day of year.
 *
 * @param[in]  year year number ussually 4 digits
 * @param[in]  month month number in [1,12]
 * @param[in]  day_of_month day of the month in [1,31]
 * @param[out] day_of_year integer day number of the year in [1,366]
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_ymd_to_day_of_year(int year, int month, int day_of_month,
        int * const day_of_year);


/**
 * The procedure "julian_to_date" does the reverse operation of the procedure
 * "make_julian_day" i.e. computes the month number and the respective day of
 * month from the information on year and integer day number.
 *
 * Source : M. Albuison
 *
 * @param[in]  year year number (4 digits)
 * @param[in]  day_of_year integer day number of the year in [1,366]
 * @param[out] month_number month number in [1,12]
 * @param[out] day_of_month day of the month in [1,31]
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_day_of_year_to_ymd(int year, int day_of_year, int *month_number,
        int *day_of_month);


/**
 * The procedure "nbdays_month" gives the number of days in a month, useful for
 * monthly calculations. Returns 0 if OK, 1 otherwise.
 *
 * @param[in]  year year number (4 digits)
 * @param[in]  month month number in [1,12]
 * @param[out] number_days_of_month number of days in a month in [28,31]
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_nbdays_month(int year, int month, int *number_days_of_month);


/**
 * The procedure "number_to_name_month" converts the month number into the
 * corresponding month name.
 *
 * @param[in]  month month number in [1,12]
 * @param[out] month_name French name of month (3 characters only, jan..dec)
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_number_to_name_month(int month_number, char *month_name);


/**
 * The procedure "Day_Angle" expresses the integer day number as an angle (in radians)
 * from 12:00 hours on the day 31st December. A year length of 365.2422 days is used.
 * Returns 0 if OK, 1 otherwise.
 *
 * @param[in]  day_of_year the day number within the year in [1,366]
 * @param[out] day_angle
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_day_angle(int day_of_year, double *day_angle);


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
 * @param[in]  year the year number usualy 4 digits
 * @param[in]  day_of_year the number of the day within year
 * @param[in]  lambda the longitude of the cite in radians
 * @return     declination of the sun in radians
 **/
int sg1_declination_sun(int year, int day_of_year, double lambda,
        double *delta);


/**
 * The procedure "declination_sun_month" computes the noon solar declination
 * (in radians) in solar time, with a simplified form for estimating monthly
 * mean global solar radiation.
 *
 * Source : Gruter (Ed.) (1984)
 *
 * @param[in]  month_number month number in [1,12]
 * @param[out] delta_month solar declination angle (in radians)
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_declination_sun_month_avg(int month_number, double *delta_month);


/**
 * The procedure "declination_sun_month" computes the noon solar declination
 * (in radians) in solar time, with a simplified form for estimating monthly
 * mean maximum global solar radiation.
 *
 * Source : Gruter (Ed.) (1984)
 *
 * @param[in]  month_number month number in [1,12]
 * @param[out] delta_month solar declination angle (in radians)
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_declination_sun_month_max(int month_number, double *delta_month);


/**
 * The procedure "solar_hour_angle" supplies the solar hour angle (in radians). By
 * convention the hour angle is negative before noon and positive after noon.
 *
 * @param[in]  t: solar time i.e. TST (0..24 decimal hours)
 * @param[out] omega: solar hour angle in radians.
 * @return 0 if everything is OK, otherwise return 1.
 **/
int sg1_solar_hour_angle(double t, double *omega);


/**
 * The procedure "omega_to_TST" does the reverse operation of the procedure
 * "solar_hour_angle" i.e. computes the solar time (in decimal hours) from the
 * solar hour angle (in radians).
 *
 * @param[in]  omega solar hour angle (in radians)
 * @param[out] t solar time i.e. TST decimal hours in [0,24]
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_omega_to_TST(double omega, double *t);

/**
 * Convert geograpique latitude to geocentrique one
 *
 * @param[in]  phi_g latitude geographic in radian
 * @return     latitude geocentric in radian
 **/
double sg1_geogr_to_geoce(double phi_g);


/**
 * The procedure "solar_hour_angle_h" supplies an average value of the solar
 * hour angle (in radians) for a whole solar hour, taking into account only the
 * portion of the solar hour with the sun standing above the horizon.
 *
 * @param[in]  phi_g latitude geographic in radians, positive to North
 * @param[in]  delta solar declination angle in radians
 * @param[in]  t solar time i.e. TST decimal hours in [0,24]
 * @param[out] omega solar hour angle in radians
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_solar_hour_angle_h(double phi_g, double delta, double t, double *omega);


/*********************************
 * SUNRISE, SUNSET AND DAYLENGTH *
 *********************************/


/**
 * The procedure "sunrise_hour_angle" supplies the sunrise and sunset hour
 * angles (in radians). Due to the dimension of the solar disk and the effect
 * of the atmospheric refraction, the edge of the solar disk will just appear
 * (disappear) at the horizon at sunrise (at sunset) when the calculated
 * astronomical elevation is 50'.
 *
 * @param[in]  phi_g latitude geographic in radians, positive to North
 * @param[in]  delta solar declination angle in radians
 * @param[in]  gamma_riset solar elevation near sunrise/sunset:
 *               - set to  0.0 for astronomical sunrise/sunset
 *               - set to -1.0 for refraction corrected sunrise/sunset.
 * @param[out] omega_sr sunrise solar hour angle in radians
 * @param[out] omega_ss sunset solar hour angle in radians
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_sunrise_hour_angle(double phi_g, double delta, double gamma_riset,
        double *omega_sr, double *omega_ss);

/**
 * Compute the astronomical sunset
 *
 * @param[in]  phi geocentric phi in radians (latitude)
 * @param[in]  delta sun declination in radians
 * @return omega at sunset
 **/
double sg1_sunset(double phi, double delta);


/**
 * The procedure "timerise_daylength" supplies the times of astronomical
 * sunrise and sunset, and the astronomical daylength, all in TST decimal
 * hours.
 *
 * @param[in]  omega_sr sunrise hour angle in radians
 * @param[in]  omega_ss sunset hour angle in radians
 * @param[out] t_sr time of astronomical sunrise in decimal hours
 * @param[out] t_ss time of astronomical sunset in decimal hours
 * @param[out] S0 astronomical daylength in decimal hours
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_timerise_daylength(double omega_sr, double omega_ss, double *t_sr,
        double *t_ss, double *S0);

/****************************
 * CHANGING THE TIME SYSTEM *
 ****************************/


/**
 * The procedure "LMT_to_TST computes the difference (in decimal hours) between
 * the TST (true solar time) and the LMT (local mean time or clock time)
 * systems at solar noon. Two stages:
 *   - the first stage calculates the equation of time, ET, wich allows for
 *     perturbations in the rotational and angular orbital speed of the Earth.
 *   - the second stage handles the difference between the longitude of the site
 *     under consideration and the reference time zone longitude for the site. A
 *     summer time correction must be added for some countries.
 *
 * Source : Gruter (ed.) (1984)
 *
 * @param[in]  day_angle day angle in radians
 * @param[in]  lambda longitude of the site in radians, positive to East
 * @param[in]  lambda_ref reference longitude of the time zone in radians
 * @param[in]  summer_corr correction for summer time in integral hours
 * @param[out] dt Offset between local mean time (LMT) and local apparent time
 *             (TST) in decimal hours
 * @return     Returns 0 if OK, 1 otherwise.
 */
int sg1_LMT_to_TST(double day_angle, double lambda, double lambda_ref,
        int summer_corr, double *dt);


/**
 * The procedure "UT_to_TST computes the conversion of the UT (Universal time)
 * into the TST (true solar time) systems at solar noon (in decimal hours).
 * First, the equation of time, ET, is computed (in decimal hours), wich allows
 * for perturbations in the rotational and angular orbital speed of the Earth.
 *
 * @param[in]  UT hours Universal Time
 * @param[in]  day_angle day angle in radians
 * @param[in]  lambda longitude of the site in radians, positive to East
 * @param[out] TST local apparent time or solar time or true solar time (TST)
 *             in decimal hours
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_UT_to_TST(double UT, double day_angle, double lambda, double *TST);

/**********************************
 * POSITION OF THE SUN IN THE SKY *
 **********************************/


/**
 * The procedure "elevation_zenith_sun" computes the solar elevation (or
 * altitude) angle and the solar zenithal (or incidence) angle. These two
 * angles are complementary.
 *
 * @param[in]  phi_g latitude geographic of site in radians, positive to North
 * @param[in]  delta solar declination angle in radians
 * @param[in]  omega solar hour angle in radians
 * @param[out] gamma solar altitude angle in radians
 * @param[out] theta solar zenithal angle in radians
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_elevation_zenith_sun(double phi_g, double delta, double omega,
        double *gamma, double *theta);


/**
 * The procedure "azimuth_sun" computes the solar azimuth angle in the Northern
 * hemisphere. The azimuth angle has a positive value when the sun is to the
 * west of South, i.e. during the afternoon in solar time. For the Southern
 * hemisphere, the azimuth angle is measured from North. Returns 0 if OK, 1
 * otherwise.
 *
 * @param[in]  phi_g latitude geographic of site in radians, positive to North
 * @param[in]  delta solar declination angle in radians
 * @param[in]  omega solar hour angle in radians
 * @param[in]  gamma solar altitude angle in radians
 * @param[out] alpha : solar azimuthal angle in radians
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_azimuth_sun(double phi_g, double delta, double omega, double gamma,
        double *alpha);

/********************************
 * EXTRATERRESTRIAL IRRADIATION *
 ********************************/


/**
 * The procedure "corr_distance" computes the correction for the variation of
 * sun-earth distance from its mean value (also known as eccentricity). It is a
 * fucntion of time, but a single (average) value per day is enough for
 * practical calculations.
 *
 * Source : Gruter (ed.) (1984)
 *
 * @param[in]  day_angle day angle in radians
 * @param[out] eccentricity correction for Earth orbit eccentricity
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_corr_distance(double day_angle, double *eccentricity);


/**
 * The procedure "G0_normal" delivers the extraterrestrial solar irradiance
 * normal to beam for day j.
 *
 * @param[in]  IOj extraterrestrial solar irradiance normal to beam for day j;
 *                 I0j=I0*fj
 * @param[in]  theta solar incidence angle or solar zenithal angle
 * @param[out] G0 extraterrestrial global solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_G0_normal(double I0j, double theta, double *G0);


/**
 * The procedure "G0_general" delivers the extraterrestrial solar irradiation
 * incident on an horizontal surface in the general case in Wh/m2.
 *
 * @param[in]  phi_g latitude geographic of site in radians, positive to North
 * @param[in]  eccentricity correction for Earth orbit eccentricity
 * @param[in]  delta solar declination angle in radians
 * @param[in]  omega1 solar hour angle at beginning of the period in radians
 * @param[in]  omega2 solar hour angle at end of the period in radians
 * @param[out] G0_12 : extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_G0_general(double phi_g, double eccentricity, double delta, double omega1,
        double omega2, double *G0_12);


/**
 * The procedure "G0_day" delivers the extraterrestrial solar irradiation
 * incident on an horizontal surface in case of daily values (in Wh/m2), i.e.
 * omega1 = omega_sr = -omega_ss  et omega2 = omega_ss.
 *
 * REMARK: It is a special case of G0_general with the sunrise and sunset
 * angles as integration limits.
 *
 * @param[in]  phi_g latitude geographic of site in radians, positive to North
 * @param[in]  eccentricity correction for Earth orbit eccentricity
 * @param[in]  delta solar declination angle in radians
 * @param[out] G0d : daily extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_G0_day(double phi_g, double eccentricity, double delta, double *G0d);


/**
 * The procedure "G0_hours_profile" delivers the extraterrestrial solar
 * irradiation incident on an horizontal surface in case of hourly values, for
 * the 24 integral hours in a given day in Wh/m2, i.e. |omega1-omega2| =
 * Pi/12.
 *
 * @param[in]  phi_g latitude geographic of site in radians, positive to North
 * @param[in]  eccentricity correction for Earth orbit eccentricity
 * @param[in]  delta solar declination angle in radians
 * @param[out] G0h[24] 24 hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_G0_hours_profile(double phi_g, double eccentricity, double delta,
        double *G0h);


/**
 * The procedure "G0_hour" delivers the extraterrestrial solar irradiation
 * incident on an horizontal surface for a specific hour in a given day in
 * Wh/m2, i.e. |omega1-omega2| = Pi/12. t is taken as the mid hour for
 * computation of the hourly value of extraterrestrial solar irradiation.
 *
 * @param[in]  phi_g latitude geographic of site in radians, positive to North
 * @param[in]  eccentricity correction for Earth orbit eccentricity
 * @param[in]  delta solar declination angle in radians
 * @param[in]  t solar time i.e. LAT (0..24 decimal hours)
 * @param[out] G0h hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_G0_hour(double phi_g, double eccentricity, double delta, double t,
        double *G0h);


/***********************************************
 * MONTHLY AVERAGES OF SOLAR INPUTS PARAMETERS *
 ***********************************************/


/**
 * The procedure "monthly_averages" computes directly the monthly average
 * values of solar parameters : day angle (in radians), eccentricity,
 * declination (in radians), sunset hour angle (in radians), daylength (in
 * decimal hours), daily extraterrestrial irradiation (in Wh/m2) and the 24
 * hourly extraterrestrial solar irradiation (in Wh/m2).
 *
 * @param[in]  month_number month number in [1,12]
 * @param[in]  year_number year number ussually 4 digits
 * @param[in]  phi_g latitude of site in radians, positive to North
 * @param[in]  lambda longitude of site in radians, positive to East
 * @param[in]  gamma_riset solar elevation near sunrise/sunset:
 *               - set to  0.0 for astronomical sunrise/sunset
 *               - set to -1.0 for refraction corrected sunrise/sunset.
 * @param[out] day_angle_m day angle in radians
 * @param[out] delta_m solar declination angle in radians
 * @param[out] omega_ss_m sunset hour angle in radians
 * @param[out] S0_m astronomical daylength in decimal hours
 * @param[out] eccentricity_m eccentricity
 * @param[out] G0d_m daily extraterrestrial irradiation in Wh/m2
 * @param[out] G0h_m[24] 24 hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_monthly_averages(int month_number, int year_number, double phi_g,
        double lambda, double gamma_riset, double *day_angle_m, double *delta_m,
        double *omega_ss_m, double *S0_m, double *eccentricity_m, double *G0d_m,
        double *G0h_m);


/****************************************************
 * YEARLY AVERAGES OF MONTHLY SOLAR PARAMETERS      *
 * (LONG TERM MEANS OF MONTHLY MEANS OF DAILY SUMS) *
 ****************************************************/


/**
 * The procedure "yearly_averages" computes directly the yearly average over a
 * defined period of years of monthly average values of solar parameters : day
 * angle in radians, eccentricity, declination in radians, sunset hour
 * angle in radians, daylength in decimal hours, daily extraterrestrial
 * irradiation in Wh/m2 and 24 hourly extraterrestrial solar irradiation
 * in Wh/m2.
 *
 * @param[in]  month_number month number in [1,12]
 * @param[in]  year_start starting year of the considered period ussually 4 digits
 * @param[in]  year_end ending year of the considered period ussually 4 digits
 * @param[in]  phi_g latitude of site in radians, positive to North
 * @param[in]  lambda longitude of site in radians, positive to East
 * @param[in]  gamma_riset solar elevation near sunrise/sunset:
 *               - set to  0.0 for astronomical sunrise/sunset
 *               - set to -1.0 for refraction corrected sunrise/sunset.
 * @param[out] day_angle_y day angle in radians
 * @param[out] delta_y solar declination angle in radians
 * @param[out] omega_ss_y sunset hour angle in radians
 * @param[out] S0_y astronomical daylength in decimal hours
 * @param[out] eccentricity_y eccentricity
 * @param[out] G0d_y daily extraterrestrial irradiation in Wh/m2
 * @param[out] G0h_y[24] 24 hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_yearly_averages(int month_number, int year_start, int year_end,
        double phi_g, double lambda, double gamma_riset, double *day_angle_y,
        double *delta_y, double *omega_ss_y, double *S0_y,
        double *eccentricity_y, double *G0d_y, double *G0h_y);


/*********************************************
 * SOLAR INPUTS PARAMETERS FOR A CERTAIN DAY *
 *********************************************/


/**
 * The procedure "solar_parameters_day" computes the solar geometry related
 * values for a certain day : day angle (in radians), eccentricity,
 * declination (in radians), sunset hour angle (in radians), daylength (in
 * decimal hours), daily extraterrestrial irradiation (in Wh/m2) and the 24
 * hourly extraterrestrial solar irradiation (in Wh/m2).
 *
 * REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar.
 *
 * @param[in]  day_of_month day of the month in [1,31]
 * @param[in]  month_number month number in [1,12]
 * @param[in]  year_number year number  ussually 4 digits
 * @param[in]  phi_g latitude of site in radians, positive to North
 * @param[in]  lambda longitude of site in radians, positive to East
 * @param[in]  gamma_riset solar elevation near sunrise/sunset:
 *               - set to  0.0 for astronomical sunrise/sunset
 *               - set to -1.0 for refraction corrected sunrise/sunset.
 * @param[out] day_angle day angle in radians
 * @param[out] delta solar declination angle in radians
 * @param[out] omega_ss sunset hour angle in radians
 * @param[out] S0 astronomical daylength in decimal hours
 * @param[out] eccentricity eccentricity
 * @param[out] G0d daily extraterrestrial irradiation in Wh/m2
 * @param[out] G0h[24] 24 hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_solar_parameters_day(int day_of_month, int month_number, int year_number,
        double phi_g, double lambda, double gamma_riset, double *day_angle,
        double *delta, double *omega_ss, double *S0, double *eccentricity,
        double *G0d, double *G0h);


/******************************************************************
 * SOLAR INPUTS PARAMETERS FOR MONTHLY AVERAGE IRRADIATION MODELS *
 ******************************************************************/


/**
 * The procedure "solar_parameters_avg" computes the solar geometry related
 * values for monthly average irradiation models : day angle (in radians),
 * eccentricity, declination (in radians), sunset hour angle (in radians),
 * daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2)
 * and the 24 hourly extraterrestrial solar irradiation (in Wh/m2).
 *
 * REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar.
 *
 * @param[in]  month_number month number in [1,12]
 * @param[in]  phi_g latitude of site in radians, positive to North
 * @param[in]  gamma_riset solar elevation near sunrise/sunset:
 *               - set to  0.0 for astronomical sunrise/sunset
 *               - set to -1.0 for refraction corrected sunrise/sunset.
 * @param[out] day_angle_avg day angle in radians
 * @param[out] delta_avg solar declination angle in radians
 * @param[out] omega_ss_avg sunset hour angle in radians
 * @param[out] S0_avg astronomical daylength in decimal hours
 * @param[out] eccentricity_avg eccentricity
 * @param[out] G0d_avg daily extraterrestrial irradiation in Wh/m2
 * @param[out] G0h_avg[24] 24 hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_solar_parameters_avg(int month_number, double phi_g, double gamma_riset,
        double *day_angle_avg, double *delta_avg, double *omega_ss_avg,
        double *S0_avg, double *eccentricity_avg, double *G0d_avg,
        double *G0h_avg);


/******************************************************************
 * SOLAR INPUTS PARAMETERS FOR MONTHLY MAXIMUM IRRADIATION MODELS *
 ******************************************************************/


/**
 * The procedure "solar_parameters_max" computes the solar geometry related
 * values for monthly average irradiation models : day angle (in radians),
 * eccentricity, declination (in radians), sunset hour angle (in radians),
 * daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2)
 * and the 24 hourly extraterrestrial solar irradiation (in Wh/m2).
 *
 * REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar.
 *
 * @param[in]  month_number month number in [1,12]
 * @param[in]  phi_g latitude of site in radians, positive to North
 * @param[in]  gamma_riset solar elevation near sunrise/sunset:
 *               - set to  0.0 for astronomical sunrise/sunset
 *               - set to -1.0 for refraction corrected sunrise/sunset.
 * @param[out] day_angle_max day angle in radians
 * @param[out] delta_max solar declination angle in radians
 * @param[out] omega_ss_max sunset hour angle in radians
 * @param[out] S0_max astronomical daylength in decimal hours
 * @param[out] eccentricity_max eccentricity
 * @param[out] G0d_max daily extraterrestrial irradiation in Wh/m2
 * @param[out] G0h_max[24] 24 hourly extraterrestrial solar irradiation in Wh/m2
 * @return     Returns 0 if OK, 1 otherwise.
 **/
int sg1_solar_parameters_max(int month_number, double phi_g, double gamma_riset,
        double *day_angle_max, double *delta_max, double *omega_ss_max,
        double *S0_max, double *eccentricity_max, double *G0d_max,
        double *G0h_max);

/**
 * TODO: documentation
 **/
int sg1_intervals_omega_tilted_plane(double phi_g, double delta, double omega_ss,
        double beta, double alpha, double *v_om, int *p_nb);

/**
 * Convert y-m-d date to julian day (number of day from -4713
 *
 * @param[in]  year
 * @param[in]  month
 * @param[in]  day_of_month
 * @return     julian day a 12h
 **/
int sg1_ymd_to_julian_day(int year, int month, int day_of_month);


/**
 * Compute year-month-day value from a given julian day
 *
 * @param[in]  jd the julian day
 * @param[out] year the year number usually 4 digits
 * @param[out] month the month number in [1,12]
 * @param[out] day_of_month the day number within the month in [1,31]
 **/
void sg1_julian_day_to_ymd(int jd, int * year, int * month, int * day_of_month);


/***************************************/
/* POSITION OF THE SUN IN THE SKY FAST */
/***************************************/
typedef struct sg1_sgf {

    double phi_g;
    double phi;
    double delta;

    double sin_phi;
    double sin_delta;
    double cos_phi;
    double cos_delta;
    double sin_phi_sin_delta;
    double cos_phi_cos_delta;

    /* For computation of the cosinus of the incident angle */

    double alpha; //< tilted plan azimuth
    double beta; //< tilted plan tilt, or inclination

    double cos_alpha;
    double sin_alpha;

    double cos_beta;
    double sin_beta;

    double A;
    double B;
    double C;
} SG1_SOLAR_GEOMETRY_FAST;

void sg1_init_solar_geometry_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double phi_g,
        double delta);
void sg1_deftilt_solar_geometry_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double alpha,
        double beta);
void sg1_elevation_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega,
        double *p_gamma);
void sg1_elevation_zenith_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega,
        double *p_gamma, double *p_theta);
void sg1_azimuth_sun_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double sin_omega,
        double gamma, double *p_alpha);
void sg1_cos_incident_angle_fast(SG1_SOLAR_GEOMETRY_FAST *p_sgf, double cos_omega,
        double sin_omega, double *p_costhetai);

#ifdef	__cplusplus
}
#endif

#ifdef  __cplusplus

/* The C++ API */

/*
 * NOTE ABOUT C++ API:
 *
 * The C++ API is intensionnaly UNSAFE, by unsafe, We mean that that API do not
 * check any input parameters to check if they are in valid range. If you want
 * safer API use the obsolete C API. Anyway, we recommand to use the C++ API.
 *
 */

#include <tuple>
#include <cmath>

namespace sg1 {

static constexpr double const PI_LOW_PRECISION = 3.141592654;
static constexpr double const I0 = 1367.0; /* solar constant in W/m2 */
static constexpr double const DAY_LENGTH = 24.0; /* average value for the length of the day in decimal hours */

static constexpr double const PI = std::acos(-1.0);
static constexpr double const PI_2 = 1.57079632679489661923;

/* convenient precalcul of omega range to accelerate irradiation computation */
class omega_range
{
	double _omega1;
	double _omega2;
	double _d_omega; /* = omega2 - omega1 */
	double _d_sin_omega; /* = sin(omega2) - sin(omega1) */

public:
	omega_range() = default;
	omega_range(omega_range const &) = default;
	omega_range & operator=(omega_range const &) = default;

	/**
	 * This function updates the period parameters.
	 * Inputs:
	 *  ths : the output struct
	 *  omega1 : hour angle at the begining of the period
	 *  omega2 : hour angle at the end of the period
	 *  omega1 must <= omega2.
	 */
	omega_range(double omega1, double omega2)
	{
		_omega1 = omega1;
		_omega2 = omega2;
		_d_omega = omega2 - omega1;
		_d_sin_omega = sin(omega2) - sin(omega1);
	}

	double omega1() const {
		return _omega1;
	}

	double omega2() const {
		return _omega2;
	}

	double d_omega() const {
		return _d_omega;
	}

	double d_sin_omega() const {
		return _d_sin_omega;
	}

};

class G0_general_daily_func
{
	double _omega_sr;
	double _omega_ss;
    double _sin_phi_x_sin_delta;
    double _cos_phi_x_cos_delta;
    double _dt_x_i0_x_eccentricity;

    double _unsafe_exec(omega_range const & range) const {
		return _dt_x_i0_x_eccentricity * (_sin_phi_x_sin_delta * range.d_omega() + _cos_phi_x_cos_delta * range.d_sin_omega());
    }

public:
    G0_general_daily_func(double phi, double delta, double eccentricity, double omega_sr, double omega_ss)
	{
    	_omega_sr = omega_sr;
    	_omega_ss = omega_ss;
		_sin_phi_x_sin_delta = std::sin(phi) * std::sin(delta);
		_cos_phi_x_cos_delta = std::cos(phi) * std::cos(delta);
		_dt_x_i0_x_eccentricity = sg1::I0 * eccentricity * sg1::DAY_LENGTH / (2.0 * sg1::PI_LOW_PRECISION);
	}

    /**
     * Compute the irradiation for the given omega range
     **/
    double operator() (omega_range const & range) const
    {
    	// FIXME: issue when omega1 or omega2 are not in [-pi,pi]
    	if (std::isnan(_omega_sr))
    		return 0.0;
    	if (std::isnan(_omega_ss))
    		return 0.0;
    	if (range.omega2() < _omega_sr)
    		return 0.0;
    	if (range.omega1() > _omega_ss)
    		return 0.0;

    	if (range.omega1() < _omega_sr or range.omega2() > _omega_ss) {
    		double omega1 = range.omega1() < _omega_sr ? _omega_sr : range.omega1();
    		double omega2 = range.omega1() > _omega_ss ? _omega_ss : range.omega2();
    		omega_range fix_range(omega1, omega2);
    		return _unsafe_exec(fix_range);
    	} else {
    		return _unsafe_exec(range);
    	}
    }

    /**
     * Compute the irradiation for the whole day
     * equivalent to G0_general_dayly_func(omega_range{omega_sr, omega_ss});
     **/
    double operator()() const
    {
    	return _unsafe_exec(omega_range{_omega_sr, _omega_ss});
    }

    double omega_sr() const {
    	return _omega_sr;
    }

    double omega_ss() const {
    	return _omega_ss;
    }

    double sin_phi_x_sin_delta() const {
    	return _sin_phi_x_sin_delta;
    }

    double cos_phi_x_cos_delta() const {
    	return _cos_phi_x_cos_delta;
    }

    double dt_x_i0_x_eccentricity() const {
    	return _dt_x_i0_x_eccentricity;
    }

};


/**
 * @param[in]  phi geocentric phi in radians (latitude)
 * @param[in]  delta sun declination in radians
 * @return     omega at sunset
 **/
double omega_sunset(double phi, double delta);

/**
 * @param[in]  phi geocentric phi in radians (latitude)
 * @param[in]  delta sun declination in radians
 * @return     gamma sun
 **/
double gamma_sun(double phi, double delta, double omega);


/**
 * The procedure "Day_Angle" expresses the integer day number as an angle (in radians)
 * from 12:00 hours on the day 31st December. A year length of 365.2422 days is used.
 *
 * @param[in]  day_of_year the day number within the year in [1,366]
 * @return Day_Angle
 **/
double day_angle(int day_of_year);


/**
 * The procedure "corr_distance" computes the correction for the variation of sun-earth
 * distance from its mean value (also known as eccentricity). It is a fucntion of time,
 * but a single (average) value per day is enough for practical calculations.
 *
 * Source : Gruter (ed.) (1984)
 *
 * @param[in]  day_angle Day_Angle in radians
 * @return     eccentricity
 **/
double corr_distance(double day_angle);


/**
 * The procedure "nbdays_month" gives the number of days in a month, useful for monthly
 * calculations.
 *
 * @param[in]  year: the year number usualy in 4 digits
 * @param[in]  month: the number of the month
 **/
int nbdays_month(int year, int month);


/**
 * The procedure "declination_sun" computes the solar declination at noon in solar time
 * (in radians). A single (average) value per day -at noon- is adequate for pratical
 * calculations. The noon declination depends on longitude, as noon occurs earlier if
 * longitude is East of Greenwi²ch, and later if it is West. The chosen algorithm uses
 * 1957 as base year; it is basically a truncated Fourier series with six harmonics.
 *
 * Sources : Bourges, B., 1985. Improvement in solar declination computation. Solar
 * Energy, 35 (4), 367-369. Carvalho, M.J. and Bourges, B., 1986. Program Eufrad 2.0 -
 * User's Guide. Project EUFRAT final scientific report, Contract EN3S-0111-F, Solar
 * Energy and Development in the European Community, pp. 12.1-12.74. Duffie, J.A. and
 * Beckman, W.A., 1980. Solar Engineering of Thermal Processes. Wiley-Interscience, New
 * York.
 *
 * @param[in]  year the year number usualy 4 digits
 * @param[in]  day_of_year the number of the day within year
 * @param[in]  lambda the lingitude of the cite in radians
 * @return     declination of the sun in radians
 **/
double declination_sun(int year, int day_of_year, double lambda);


/**
 * The procedure "solar_hour_angle" supplies the solar hour angle (in radians). By
 * convention the hour angle is negative before noon and positive after noon.
 *
 * @param[in]  t solar time i.e. TST (0..24 decimal hours)
 * @return     solar_hour_angle in radians.
 **/
double solar_hour_angle(double t);


/**
 * Supplies the solar time in hours in [0,24].
 *
 * @param[in]  omega solar_hour_angle in radians
 * @return     solar time i.e. TST (0..24 decimal hours)
 **/
double omega_to_TST(double omega);

/**
 * TODO: Dococumentation
 **/
int ymd_to_day_of_year(int year, int month, int day_of_month);

/**
 * Compute year-month-day value from a given julian day
 *
 * @param[in]  jd the julian day
 * @return     year, month, day_of_month tuple.
 *             year the year number usually 4 digits
 *             month the month number in [1,12]
 *             day_of_month the day number within the month in [1,31]
 *             Tips: use std::tie(year, month, day)
 **/
std::tuple<int,int,int> julian_date_to_ymd(int jd);

/**
 * Compute julian date from year, month, day_of_month
 *
 * @param[in]  year the year number usualy 4 digits
 * @param[in]  month the month number within [1,12]
 * @param[in]  day_of_month the day number of the month within [1,31]
 * @return     julian date coresponding to the year, month, day_of_month provide
 **/
int ymd_to_julian_date(int year, int month, int day_of_month);

/**
 * The procedure "day_of_year_to_ymd" does the reverse operation of the procedure
 * "ymd_to_day_of_year" i.e. computes the month number and the respective day of
 * month from the information on year and integer day number.
 *
 * @param[in]  year year number (4 digits)
 * @param[in]  day_of_year integer day number of the year in [1,366]
 * @return     month, day_of_month tuple.
 *             month the month number in [1,12]
 *             day_of_month the day number within the month in [1,31]
 *             Tips: use std::tie(month, day_of_month)
 **/
std::tuple<int, int> day_of_year_to_ymd(int year, int day_of_year);

/**
 * TODO: Documentation
 **/
double geogr_to_geoce(double phi_g);

/**
 * TODO: Documentation
 **/
double azimuth_sun(double phi, double delta, double omega, double gamma);

/**
 * TODO: Documentation
 **/
double G0_general(double phi_g, double eccentricity, double delta,
		double omega1, double omega2);


/**
 * The procedure "LMT_to_TST computes the difference (in decimal hours) between
 * the TST (true solar time) and the LMT (local mean time or clock time)
 * systems at solar noon. Two stages:
 *   - the first stage calculates the equation of time, ET, wich allows for
 *     perturbations in the rotational and angular orbital speed of the Earth.
 *   - the second stage handles the difference between the longitude of the site
 *     under consideration and the reference time zone longitude for the site. A
 *     summer time correction must be added for some countries.
 *
 * Source : Gruter (ed.) (1984)
 *
 * @param[in]  day_angle day angle in radians
 * @param[in]  lambda longitude of the site in radians, positive to East
 * @param[in]  lambda_ref reference longitude of the time zone in radians
 * @param[in]  summer_corr correction for summer time in integral hours
 * @param[out] dt Offset between local mean time (LMT) and local apparent time
 *             (TST) in decimal hours
 * @return     Returns 0 if OK, 1 otherwise.
 */
double lmt_to_tst(double day_angle, double lambda, double lambda_ref, int summer_corr);


/**
 * The procedure "ut_to_tst computes the conversion of the UT (Universal time)
 * into the TST (true solar time) systems at solar noon (in decimal hours).
 * First, the equation of time, ET, is computed (in decimal hours), wich allows
 * for perturbations in the rotational and angular orbital speed of the Earth.
 *
 * @param[in]  UT hours Universal Time
 * @param[in]  day_angle day angle in radians
 * @param[in]  lambda longitude of the site in radians, positive to East
 * @return     TST local apparent time or solar time or true solar time (TST)
 *             in decimal hours
 **/
double ut_to_tst(double ut, double day_angle, double lambda);

} // namespace sg1

#endif // C++ API

#endif // __H_solar_geometry
