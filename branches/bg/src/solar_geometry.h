/*-------------------------------------------------------------------------*/
/*              GEOMETRY OF SOLAR BEAM, A NEAR POINTSOURCE                 */
/*-------------------------------------------------------------------------*/
/*                        ECOLE DES MINES DE PARIS                         */
/*        CENTRE D'ENERGETIQUE - GROUPE TELEDETECTION & MODELISATION       */
/*                       Rue Claude Daunesse, BP 207                       */
/*                   06904 Sophia Antipolis cedex, FRANCE                  */
/*          Tel (+33) 04 93 95 74 49     Fax (+33) 04 93 95 75 35          */
/*                  E-mail : lucien.wald@mines-paristech.fr                */
/*-------------------------------------------------------------------------*/
/*   L. Wald - O. Bauer - February 1997                                    */
/*   modified 8 July 2004 L. Wald for geocentric - geographic lat          */
/*-------------------------------------------------------------------------*/

#ifndef __H_solar_geometry
#define __H_solar_geometry

#ifdef _MINGW_
#include <windows.h>
#endif

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef PUBLIC
#undef PUBLIC
#endif

#ifdef INIT
#undef INIT
#endif

#ifdef __C_solar_geometry
#define INIT
#define PRIVATE static
#endif

#ifdef _VISUAL_
#define EXPORT _declspec (dllexport)
#else
#define EXPORT
#endif

#define I0  1367.0		/* solar constant in W/m2 */
#define DAY_LENGTH  24.0		/* average value for the length of the day in decimal hours */

/*
 * NOTATIONS
 * ----------
 *
 * phi_g : geographic latitude of the site, positive to North
 * phi : geocentric latitude of the site, positive to North
 * lambda : longitude of the site, positive to East
 * delta : solar declination angle
 * omega : solar hour angle
 * gamma : solar altitude (or elevation) angle
 * theta : solar incidence (or zenithal) angle (ESRA --> zeta)
 * alpha : solar azimuthal angle (or psi)
 * t : solar time = true solar time (TST) = local apparent time (LAT)
 * LAT : local apparent time or solar time or true solar time (TST)
 *  --> this system of time offers the advantage of symmetry of the solar geometry about the north-south line
 * LMT : local mean time or clock time
 * UT : Universal Time, is GMT measured from Greenwich mean midnight
 * omega_sr : sunrise hour angle
 * omega_ss : sunset hour angle
 * t_sr : time of astronomical sunrise
 * t_ss : time of astronomical sunset
 * omega1 : solar hour angle at beginning of the time period
 * omega2 : solar hour angle at end of the time period
 * S0 : astronomical daylength or astronomical sunshine duration
 * I0 : solar constant = annual mean value of extraterrestrial direct solar irradiance G0 (1367.0 W/m2)
 * G0 : extraterrestrial global solar irradiation (on an horizontal plane) =B0
 * G0h : hourly extraterrestrial solar irradiation (on an horizontal plane)
 * G0d : daily extraterrestrial solar irradiation (on an horizontal plane)
 *
 * NB : All angles are computed in radians as a standard. The basic trigonometric
 * calculations on the position of the sun are carried out in LAT.
 */

/*
 * 0. BASIC PARAMETERS
 * --------------------
 */

/**
 * The procedure converts a day given in year-month-day into a
 * day of year.
 * Source :
 * @param day_of_month day of the month (1..31)
 * @param month_number month number (1..12)
 * @param year_number year number (4 digits)
 * @return day of year (1..366)
 */
EXPORT int ymd_to_day_of_year(int year, int month, int day_of_month);

/**
 * The procedure does the reverse operation of the procedure
 * "make_julian_day" i.e. computes the month number and the respective day of month from
 * the information on year and integer day number. Returns 0 if OK, 1 otherwise.
 * Source : Michel Albuisson in /u2/tm/src/srcgeo/julian_lib/
 * @param year_number year number (4 digits)
 * @param day_of_day integer day number in year (1..366)
 * @param day_of_month the resulting day of the month (1..31)
 * @param month_number the resulting month number (1..12)
 */
EXPORT void day_of_year_to_ymd(int year_number, int day_of_year,
		int * day_of_month, int * month_number);

/**
 * The procedure "nbdays_month" gives the number of days in a month, useful for monthly
 * calculations. Returns 0 if OK, 1 otherwise.
 * Source :
 * @param year_number year number (4 digits)
 * @param month_number month number (1..12)
 * @return number_days_month number of days in a month
 */
EXPORT int nbdays_month(int year_number, int month_number);

/**
 * The procedure "get_month_name" converts the month number into the corresponding
 * month name. Returns 0 if OK, 1 otherwise.
 * Source :
 * @param month_number month number (1..12)
 * @param month_name name of month (3 characters only, jan..dec)
 * @return month_name a pointer to thename of the month abbreviated with 3 characters (jan..dec)
 */
EXPORT char const * const get_month_name(int month);

/**
 * The procedure "get_day_angle" expresses the integer day number as an angle (in radians)
 * from 12:00 hours on the day 31st December. A year length of 365.2422 days is used.
 * Returns 0 if OK, 1 otherwise.
 * Source :
 * @param julian_day integer day number or julian day (1..366)
 * @return day angle (in radians)
 */
EXPORT double get_day_angle(int day_of_year);

/**
 * The procedure "declination_sun" computes the solar declination at noon in solar time
 * (in radians). A single (average) value per day -at noon- is adequate for pratical
 * calculations. The noon declination depends on longitude, as noon occurs earlier if
 * longitude is East of Greenwitch, and later if it is West. The chosen algorithm uses
 * 1957 as base year; it is basically a truncated Fourier series with six harmonics.
 * Returns 0 if OK, 1 otherwise.
 * Sources : Bourges, B., 1985. Improvement in solar declination computation. Solar
 * Energy, 35 (4), 367-369. Carvalho, M.J. and Bourges, B., 1986. Program Eufrad 2.0 -
 * User's Guide. Project EUFRAT final scientific report, Contract EN3S-0111-F, Solar
 * Energy and Development in the European Community, pp. 12.1-12.74. Duffie, J.A. and
 * Beckman, W.A., 1980. Solar Engineering of Thermal Processes. Wiley-Interscience, New
 * York.
 * @param year_number year number (4 digits)
 * @param julian_day integer day number in year (1..366)
 * @param lambda longitude (in radians, positive to East)
 * @return solar declination angle at noon (in radians)
 */
EXPORT double declination_sun(int year, int day_of_year, double lambda);

/**
 * The procedure "declination_sun_month" computes the noon solar declination (in
 * radians) in solar time, with a simplified form, in two cases: type_use=0 : for
 * estimating monthly mean global solar radiation type_use=1 : for estimating monthly
 * mean maximum global solar radiation The integer day number to be selected in each case
 * for the computations is given by two tables. Returns 0 if OK, 1 otherwise.
 * Source : Gruter (Ed.) (1984)
 * @param month_number month number (1..12)
 * @return solar declination angle (in radians)
 */
EXPORT double declination_sun_month(int month, int type_use);

/**
 * The procedure "solar_hour_angle" supplies the solar hour angle (in radians). By
 * convention the hour angle is negative before noon and positive after noon Returns 0 if
 * OK, 1 otherwise.
 * Source :
 * @param t solar time i.e. LAT (0..24 decimal hours)
 * @return solar hour angle (in radians)
 */
EXPORT double solar_hour_angle(double t);

/**
 * The procedure "omega_to_LAT" does the reverse operation of the procedure
 * "solar_hour_angle" i.e. computes the solar time (in decimal hours) from the solar
 * hour angle (in radians). Returns 0 if OK, 1 otherwise.
 * Source :
 * @param omega solar hour angle (in radians)
 * @return solar time i.e. LAT (0..24 decimal hours)
 */
EXPORT double omega_to_LAT(double omega);

/**
 * Convert geographic to geocentric
 * @param phi_g latitude geographic in radian
 * @return latitude geocentric in radian
 */
EXPORT double geogr_to_geoce(double phi_g);

/**
 * The procedure "solar_hour_angle_h" supplies an average value of the solar hour angle
 * (in radians) for a whole solar hour, taking into account only the portion of the solar
 * hour with the sun standing above the horizon. Returns 0 if OK, 1 otherwise.
 * Source :
 * @param phi_g latitude (in radians, positive to North)
 * @param delta solar declination angle (in radians)
 * @param t solar time i.e. LAT (0..24 decimal hours)
 * @return solar hour angle (in radians)
 */
EXPORT double solar_hour_angle_h(double phi_g, double delta, double t);

/*********************************/
/* SUNRISE, SUNSET AND DAYLENGTH */
/*********************************/

/**
 * The procedure "sunrise_hour_angle" supplies the sunrise and sunset hour angles (in
 * radians). Due to the dimension of the solar disk and the effect of the atmospheric
 * refraction, the edge of the solar disk will just appear (disappear) at the horizon at
 * sunrise (at sunset) when the calculated astronomical elevation is 50'. Returns 0 if
 * OK, 1 otherwise.
 * Source :
 *
 * @param phi_g latitude of site (in radians, positive to North)
 * @param delta solar declination angle (in radians)
 * @param gamma_riset solar elevation near sunrise/sunset: -
 * set to 0.0 for astronomical sunrise/sunset - set to -1.0 for refraction corrected
 * sunrise/sunset.
 * @param omega_sr the resulting sunrise solar hour angle (in radians)
 * @param omega_ss the resulting sunset solar hour angle (in radians)
 */
EXPORT void sunrise_hour_angle(double phi_g, double delta, double gamma_riset,
		double *omega_sr, double *omega_ss);

/**
 * The procedure "timerise_daylength" supplies the times of astronomical sunrise and
 * sunset, and the astronomical daylength, all in LAT decimal hours. Returns 0 if OK, 1
 * otherwise.
 * Source :
 * @param omega_sr sunrise hour angle (in radians)
 * @param omega_ss sunset hour angle (in radians)
 * @param t_sr the resulting time of astronomical sunrise (in decimal hours)
 * @param t_ss the resulting time of astronomical sunset (in decimal hours) S0 : astronomical daylength (in decimal hours)
 */
EXPORT void timerise_daylength(double omega_sr, double omega_ss, double *t_sr,
		double *t_ss, double *S0);

/*
 * 2. CHANGING THE TIME SYSTEM
 * ----------------------------
 */

/** The procedure "LMT_to_LAT computes the difference (in decimal hours) between the LAT
 * (local apparent time) and the LMT (local mean time or clock time) systems at solar
 * noon. Two stages: - the first stage calculates the equation of time, ET, wich allows
 * for perturbations in the rotational and angular orbital speed of the Earth. - the
 * second stage handles the difference between the longitude of the site under
 * consideration and the reference time zone longitude for the site. A summer time
 * correction must be added for some countries. Returns 0 if OK, 1 otherwise.
 * Source : Gruter (ed.) (1984)
 * @param day_angle day angle (in radians)
 * @param lambda longitude of the site (in radians, positive to East)
 * @param lambda_ref reference longitude of the time zone (in radians)
 * @param summer_corr correction for summer time (integer hours)
 * @return  Offset between local mean time (LMT) and local apparent time (LAT) (in decimal hours)
 */
EXPORT double LMT_to_LAT(double day_angle, double lambda, double lambda_ref,
		int summer_corr);

/**
 * The procedure "UT_to_LAT computes the conversion of the UT (Universal time) into the
 * LAT (local apparent time) systems at solar noon (in decimal hours). First, the
 * equation of time, ET, is computed (in decimal hours), wich allows for perturbations in
 * the rotational and angular orbital speed of the Earth.
 * Source :
 * @param UT Universal Time (in decimal hours)
 * @param day_angle day angle (in radians)
 * @param lambda longitude of the site (in radians, positive to East)
 * @return local apparent time or solar time or true solar time (TST) (in decimal hours)
 */
EXPORT double UT_to_LAT(double UT, double day_angle, double lambda);

/*
 * 3. POSITION OF THE SUN IN THE SKY
 * ----------------------------------
 */

/**
 * The procedure "elevation_zenith_sun" computes the solar elevation (or altitude) angle
 * and the solar zenithal (or incidence) angle. These two angles are complementary.
 * Source :
 * @param phi_g latitude of site (in radians, positive to North)
 * @param delta solar declination angle (in radians)
 * @param omega solar hour angle (in radians)
 * @param gamma the resulting solar altitude angle (in radians)
 * @param theta the resulting solar zenithal angle (in radians)
 */
EXPORT void elevation_zenith_sun(double phi_g, double delta, double omega,
		double *gamma, double *theta);

/**
 * The procedure "azimuth_sun" computes the solar azimuth angle in the Northern
 * hemisphere. The azimuth angle has a positive value when the sun is to the west of
 * South, i.e. during the afternoon in solar time. For the Southern hemisphere, the
 * azimuth angle is measured from North.
 * Source :
 * @param phi_g latitude of site (in radians, positive to North)
 * @parma delta solar declination angle (in radians)
 * @parma omega solar hour angle (in radians)
 * @param gamma solar altitude angle (in radians)
 * @return solar azimuthal angle (in radians)
 */
EXPORT double azimuth_sun(double phi_g, double delta, double omega,
		double gamma);

/*
 * 4. EXTRATERRESTRIAL IRRADIATION
 * --------------------------------
 */

/**
 * The procedure "corr_distance" computes the correction for the variation of sun-earth
 * distance from its mean value (also known as eccentricity). It is a fucntion of time,
 * but a single (average) value per day is enough for practical calculations.
 * Source : Gruter (ed.) (1984)
 * @param day_angle day angle (in radians)
 * @return correction for Earth orbit eccentricity
 */
EXPORT double corr_distance(double day_angle);

/**
 * The procedure "G0_normal" delivers the extraterrestrial solar irradiance normal to
 * beam for day j.
 * Source :
 * @param IOj extraterrestrial solar irradiance normal to beam for day j; I0j=I0*fj
 * @param theta solar incidence angle or solar zenithal angle
 * @return extraterrestrial global solar irradiation G0 (in Wh/m2)
 */
EXPORT double G0_normal(double I0j, double theta);

/**
 * The procedure "G0_general" delivers the extraterrestrial solar irradiation incident
 * on an horizontal surface in the general case (in Wh/m2).
 * Source :
 * @param phi_g latitude of site (in radians, positive to North)
 * @param eccentricity correction for Earth orbit eccentricity
 * @param delta solar declination angle (in radians)
 * @param omega1 : solar hour angle at beginning of the period (in radians)
 * @param omega2 : solar hour angle at end of the period (in radians)
 * @return extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT double G0_general(double phi_g, double eccentricity, double delta,
		double omega1, double omega2);

/**
 * The procedure "G0_day" delivers the extraterrestrial solar irradiation incident on an
 * horizontal surface in case of daily values (in Wh/m2), i.e. omega1 = omega_sr =
 * -omega_ss et omega2 = omega_ss. Returns 0 if OK, 1 otherwise. REMARK: It is a special
 * case of G0_general with the sunrise and sunset angles as integration limits.
 * Source :
 * @param phi_g latitude of site (in radians, positive to North)
 * @param eccentricity correction for Earth orbit eccentricity delta : solar declination angle (in radians)
 * @return daily extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT double G0_day(double phi_g, double eccentricity, double delta);

/**
 * The procedure "G0_hours_profile" delivers the extraterrestrial solar irradiation
 * incident on an horizontal surface in case of hourly values, for the 24 integral hours
 * in a given day (in Wh/m2), i.e. |omega1-omega2| = Pi/12.
 * Source :
 * @param phi_g latitude of site (in radians, positive to North)
 * @param eccentricity correction for Earth orbit eccentricity
 * @param delta solar declination (in radians)
 * @param the resulting G0h[1..24] : 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT void G0_hours_profile(double phi_g, double eccentricity, double delta,
		double * G0h);

/**
 * The procedure "G0_hour" delivers the extraterrestrial solar irradiation incident on
 * an horizontal surface for a specific hour in a given day (in Wh/m2), i.e.
 * |omega1-omega2| = Pi/12. t is taken as the mid hour for computation of the hourly
 * value of extraterrestrial solar irradiation. Returns 0 if OK, 1 otherwise
 * Source :
 * @param phi_g latitude of site (in radians, positive to North)
 * @param eccentricity correction for Earth orbit eccentricity
 * @param delta : solar declination (in radians)
 * @param t : solar time i.e. LAT (0..24 decimal hours)
 * @return hourly extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT double
G0_hour(double phi_g, double eccentricity, double delta, double t);

/*
 * 4.1. MONTHLY AVERAGES OF SOLAR INPUTS PARAMETERS
 * -------------------------------------------------
 */

/**
 * The procedure "monthly_averages" computes directly the monthly average
 * values of solar parameters : day angle (in radians), eccentricity,
 * declination (in radians), sunset hour angle (in radians), daylength (in
 * decimal hours), daily extraterrestrial irradiation (in Wh/m2) and the 24
 * hourly extraterrestrial solar irradiation (in Wh/m2) .
 * Source :
 * @param month month number (1..12)
 * @param year year number (4 digits)
 * @param phi_g latitude of site (in radians, positive to North)
 * @param lambda longitude of site (in radians, positive to East)
 * @param gamma_riset  : solar elevation near sunrise/sunset:
 *   - set to  0.0 for astronomical sunrise/sunset
 *   - set to -1.0 for refraction corrected sunrise/sunset.
 * @param day_angle_m the resulting monthly average of day angle (in radians)
 * @param delta_m the resulting monthly average of solar declination angle (in radians)
 * @param omega_ss_m the resulting monthly average of sunset hour angle (in radians)
 * @param S0_m the resulting monthly average of astronomical daylength (in decimal hours)
 * @param eccentricity_m the resulting monthly average of eccentricity
 * @param G0d_m the resulting monthly average of daily extraterrestrial irradiation (in Wh/m2)
 * @param G0h_m the resulting monthly average of 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT void monthly_averages(int month, int year, double phi_g, double lambda,
		double gamma_riset, double *day_angle_m, double *delta_m,
		double *omega_ss_m, double *S0_m, double *eccentricity_m,
		double *G0d_m, double *G0h_m);

/*
 * 4.2. YEARLY AVERAGES OF MONTHLY SOLAR PARAMETERS
 * (LONG TERM MEANS OF MONTHLY MEANS OF DAILY SUMS)
 * -------------------------------------------------
 */

/**
 *  The procedure "yearly_averages" computes directly the yearly average over a
 * defined period of years of monthly average values of solar parameters : day
 * angle (in radians), eccentricity, declination (in radians), sunset hour
 * angle (in radians), daylength (in decimal hours), daily extraterrestrial
 * irradiation (in Wh/m2) and 24 hourly extraterrestrial solar irradiation
 * (in Wh/m2).
 * Source :
 * @param month_number month number (1..12)
 * @param year_start starting year of the considered period (4 digits)
 * @param year_end ending year of the considered period (4 digits)
 * @param phi_g latitude of site (in radians, positive to North)
 * @param lambda longitude of site (in radians, positive to East)
 * @param gamma_riset solar elevation near sunrise/sunset:
 *   - set to  0.0 for astronomical sunrise/sunset
 *   - set to -1.0 for refraction corrected sunrise/sunset.
 * @param day_angle_y the resulting yearly average of day angle (in radians)
 * @param delta_y the resulting yearly average of solar declination angle (in radians)
 * @param omega_ss_y the resulting yearly average of sunset hour angle (in radians)
 * @param S0_y the resulting yearly average of astronomical daylength (in decimal hours)
 * @param eccentricity_y the resulting yearly average of eccentricity
 * @param G0d_y the resulting yearly average of daily extraterrestrial irradiation (in Wh/m2)
 * @param G0h_y the resulting yearly average of 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT void yearly_averages(int month, int year_start, int year_end,
		double phi_g, double lambda, double gamma_riset, double *day_angle_y,
		double *delta_y, double *omega_ss_y, double *S0_y,
		double *eccentricity_y, double *G0d_y, double *G0h_y);

/*
 * 4.3. SOLAR INPUTS PARAMETERS FOR A CERTAIN DAY
 * -----------------------------------------------
 */

/**
 * The procedure "solar_parameters_day" computes the solar geometry related
 * values for a certain day : day angle (in radians), eccentricity,
 * declination (in radians), sunset hour angle (in radians), daylength (in
 * decimal hours), daily extraterrestrial irradiation (in Wh/m2) and the 24
 * hourly extraterrestrial solar irradiation (in Wh/m2).
 * REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar.
 * Source :
 * @param day_of_month day of the month (1..31)
 * @param month_number month number (1..12)
 * @param year_number year number (4 digits)
 * @param phi_g latitude of site (in radians, positive to North)
 * @param lambda longitude of site (in radians, positive to East)
 * @param gamma_riset solar elevation near sunrise/sunset:
 *   - set to  0.0 for astronomical sunrise/sunset
 *   - set to -1.0 for refraction corrected sunrise/sunset.
 * @param day_angle day angle (in radians)
 * @param delta solar declination angle (in radians)
 * @param omega_ss sunset hour angle (in radians)
 * @param S0 astronomical daylength (in decimal hours)
 * @param eccentricity eccentricity
 * @param G0d daily extraterrestrial irradiation (in Wh/m2)
 * @param G0h 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
EXPORT void solar_parameters_day(int day_of_month, int month_number,
		int year_number, double phi_g, double lambda, double gamma_riset,
		double *day_angle, double *delta, double *omega_ss, double *S0,
		double *eccentricity, double *G0d, double *G0h);

/*
 * 4.4. SOLAR INPUTS PARAMETERS FOR MONTHLY AVERAGE IRRADIATION MODELS
 * --------------------------------------------------------------------
 */

/**
 * The procedure "solar_parameters_acg" computes the solar geometry related
 * values for monthly average irradiation models : day angle (in radians),
 * eccentricity, declination (in radians), sunset hour angle (in radians),
 * daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2)
 * and the 24 hourly extraterrestrial solar irradiation (in Wh/m2).
 * REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar.
 * Source :
 * @param month month number (1..12)
 * @param phi_g latitude of site (in radians, positive to North)
 * @param gamma_riset solar elevation near sunrise/sunset:
 *   - set to  0.0 for astronomical sunrise/sunset
 *   - set to -1.0 for refraction corrected sunrise/sunset.
 * @param day_angle_avg the resulting average of day angle (in radians) for the given month
 * @param delta_avg the resulting average of solar declination angle (in radians) for the given month
 * @param omega_ss_avg the resulting average of sunset hour angle (in radians) for the given month
 * @param S0_avg the resulting average of astronomical daylength (in decimal hours) for the given month
 * @param eccentricity_avg the resulting average of eccentricity for the given month
 * @param G0d_avg the resulting average of daily extraterrestrial irradiation (in Wh/m2) for the given month
 * @param G0h_avg the resulting average of 24 hourly extraterrestrial solar irradiation (in Wh/m2) for the given month
 */
EXPORT void solar_parameters_avg(int month, double phi_g, double gamma_riset,
		double *day_angle_avg, double *delta_avg, double *omega_ss_avg,
		double *S0_avg, double *eccentricity_avg, double *G0d_avg,
		double *G0h_avg);

/*
 * 4.5. SOLAR INPUTS PARAMETERS FOR MONTHLY MAXIMUM IRRADIATION MODELS
 * --------------------------------------------------------------------
 */

/**
 * The procedure "solar_parameters_acg" computes the solar geometry related
 * values for monthly average irradiation models : day angle (in radians),
 * eccentricity, declination (in radians), sunset hour angle (in radians),
 * daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2)
 * and the 24 hourly extraterrestrial solar irradiation (in Wh/m2).
 * REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar.
 * Source :
 * @param month_number : month number (1..12)
 * @param phi_g        : latitude of site (in radians, positive to North)
 * @param  gamma_riset  : solar elevation near sunrise/sunset:
 *   - set to  0.0 for astronomical sunrise/sunset
 *   - set to -1.0 for refraction corrected sunrise/sunset.
 * @param day_angle_max the resulting maximum of day angle (in radians) for the given month
 * @param delta_max the resulting maximum of solar declination angle (in radians) for the given month
 * @param omega_ss_max the resulting maximum of sunset hour angle (in radians) for the given month
 * @param S0_max the resulting maximum of astronomical daylength (in decimal hours) for the given month
 * @param eccentricity_max the resulting maximum of eccentricity for the given month
 * @param G0d_max the resulting maximum of daily extraterrestrial irradiation (in Wh/m2) for the given month
 * @param G0h_max the resulting maximum of 24 hourly extraterrestrial solar irradiation (in Wh/m2) for the given month
 */
EXPORT void solar_parameters_max(int month, double phi_g, double gamma_riset,
		double *day_angle_max, double *delta_max, double *omega_ss_max,
		double *S0_max, double *eccentricity_max, double *G0d_max,
		double *G0h_max);

/**
 * Compute the solar angles omega when the Sun is visible (located in front of) by the
 * tilted surface [omega is related to time TST by: tTST = 12(1+omega/pi)] an area of 4
 * values is returned: [$omega1, $omega2, $omega3, $omega4] If all values are -999, then
 * the sun is not visible by the plane at any time of the day If $omega3 and $omega4 are
 * -999, and $omega1 and $omega2 not, sun is visible during the period [$omega1;$omega2]
 * If all values are different from -999,sun is visible during the
 * periods:[$omega1;$omega2] and [$omega3;$omega4]
 * @param phi_g geographic latitude of site (in radians, positive to North)
 * @param delta solar declination angle (radian)
 * @param omega_ss : sunset hour angle (in radians)
 * @param beta tilt angle of the inclined flat plane
 * @param alpha azimuth of the inclined flat plane
 */
EXPORT void intervals_omega_tilted_plane(double phi_g, double delta,
		double omega_ss, double beta, double alpha, double *v_om, int *p_nb);

#ifdef	__cplusplus
}
#endif
#endif
