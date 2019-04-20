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
/*-------------------------------------------------------------------------*/

#ifndef __H_solar_geometry
#define __H_solar_geometry

#include "angle.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*******************************************/
/*                                         */
/* ANGLE CONVENTIONS ISO VS SOLAR GEOMETRY */
/*                                         */
/*******************************************/

/* 
 * ISO [0; 360] (independent of the hemisphere):
 * *********************************************
 *
 *                   360= NORTH = 0
 *        WEST=270                           EAST 90
 *                     SOUTH 180
 */

/*
 * Solar Geometry [-180; 180] (dependent of the hemisphere):
 * *********************************************************
 *
 *      - SG Northern Hemisphere:
 *                    180= NORTH = -180
 *        WEST = 90                           EAST = -90
 *                         SOUTH = 0
 *
 *     - SG Southern Hemisphere:
 *                         NORTH = 0
 *        WEST = 90                           EAST = -90
 *                   180 = SOUTH = -180
 */

/*************/
/* CONSTANTS */
/*************/

#define SG1_PI_LOW_PRECISION  3.141592654

#define SG1_I0  1367.0			/* solar constant in W/m2 */
#define SG1_DAY_LENGTH  24.0	/* average value for the length of the day in decimal hours */

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
/* G0  : extraterrestrial global solar irradiation (on an horizontal plane) = B0 */
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
/* BASIC PARAMETERS */
/********************/

/* Source : */
/* Inputs :
   day_of_month : day of the month (1..31)
   month_number : month number (1..12)
   year_number  : year number (4 digits) */
/* Outputs :
   julian_day : integer day number or julian day (1..366)
 */
/* The procedure "make_julian_day" converts a day given in day, month and year 
   into a julian day. Returns 0 if OK, 1 otherwise. */
int make_julian_day(int day_of_month, int month_number, int year_number,
            int *julian_day);

/* Source : MA in /u2/tm/src/srcgeo/julian_lib/ */
/* Inputs :
   year_number : year number (4 digits)
   julian_day  : integer day number or julian day (1..366)
 */
/* Outputs :
   day_of_month : day of the month (1..31)
   month_number : month number (1..12) */
/* The procedure "julian_to_date" does the reverse operation of the procedure 
   "make_julian_day" i.e. computes the month number and the respective day of 
   month from the information on year and integer day number. Returns 0 if OK,
   1 otherwise. */
int julian_to_date (int year_number, int julian_day,
		     int *day_of_month, int *month_number);

/* Source : */
/* Inputs :
   year_number  : year number (4 digits)
   month_number : month number (1..12) */
/* Outputs :
   number_days_month : number of days in a month */
/* The procedure "nbdays_month" returns the number of days in a month, useful for
   monthly calculations. */
int nbdays_month(int year_number, int month_number);

/* Source : */
/* Inputs :
   month_number : month number (1..12)
   month_name   : name of month (3 characters only, jan..dec) */
/* Outputs :
   month_name : name of the month abbreviated with 3 characters (jan..dec) */
/* The procedure "number_to_name_month" converts the month number into the 
   corresponding month name. Returns 0 if OK, 1 otherwise. */
int number_to_name_month (int month_number, char *month_name);

/* Source : */
/* Inputs :
   julian_day : integer day number or julian day (1..366) */
/* Outputs :
   day_angle : day angle (in radians) */
/* The procedure "Day_Angle" expresses the integer day number as an angle (in
   radians) from 12:00 hours on the day 31st December. A year length of 
   365.2422 days is used. */
double Day_Angle (int julian_day);

/* Sources : 
   Bourges, B., 1985. Improvement in solar declination computation. Solar 
   Energy, 35 (4), 367-369. 
   Carvalho, M.J. and Bourges, B., 1986. Program Eufrad 2.0 - User's Guide. 
   Project EUFRAT final scientific report, Contract EN3S-0111-F, Solar Energy 
   and Development in the European Community, pp. 12.1-12.74.
   Duffie, J.A. and Beckman, W.A., 1980. Solar Engineering of Thermal 
   Processes. Wiley-Interscience, New York. */
/* Inputs :
   year_number : year number (4 digits)
   julian_day  : integer day number or julian day (1..366)
   lambda      : longitude (in radians, positive to East) */
/* Outputs :
   delta : solar declination angle at noon (in radians) */
/* The procedure "declination_sun" computes the solar declination at noon in 
   solar time (in radians). A single (average) value per day -at noon- is 
   adequate for pratical calculations. The noon declination depends on 
   longitude, as noon occurs earlier if longitude is East of Greenwich, and 
   later if it is West. The chosen algorithm uses 1957 as base year; it is 
   basically a truncated Fourier series with six harmonics. */
double declination_sun(int year_number, int julian_day, double lambda);

/* Source : Gruter (Ed.) (1984) */
/* Inputs :
   month_number : month number (1..12) */
/* Outputs :
   delta_month : solar declination angle (in radians) */
/* The procedure "declination_sun_month" computes the noon solar declination 
   (in radians) in solar time, with a simplified form, in two cases:
   type_use=0 : for estimating monthly mean global solar radiation
   type_use=1 : for estimating monthly mean maximum global solar radiation
   The integer day number to be selected in each case for the computations is
   given by two tables. Returns 0 if OK, 1 otherwise. */
int declination_sun_month (int month_number, int type_use,
			    double *delta_month);

/* Source : */
/* Inputs :
   t : solar time i.e. LAT (0..24 decimal hours) */
/* Outputs :
   omega : solar hour angle (in radians) */
/* The procedure "soldar_hour_angle" supplies the solar hour angle (in radians).
   By convention the hour angle is negative before noon and positive after noon. */
double solar_hour_angle (double t);

/* Source : */
/* Inputs :
   omega : solar hour angle (in radians) */
/* Outputs :
   t : solar time i.e. LAT (0..24 decimal hours) */
/* The procedure "omega_to_LAT" does the reverse operation of the procedure 
   "solar_hour_angle" i.e. computes the solar time (in decimal hours) from the 
   solar hour angle (in radians). */
double omega_to_LAT (double omega);

/* Source : */
/* Inputs :
   phi_g : geographic latitude (in radians) */
/* Outputs :
   phi : geocentric latitude (in radians) */
/* The procedure "geogr_to_geoce" converts geographic latitude
 * into geocentric latitude. */
double geogr_to_geoce (double phi_g);

/* Source : */
/* Inputs :
   phi_g : latitude (in radians, positive to North)
   delta : solar declination angle (in radians)
   t     : solar time i.e. LAT (0..24 decimal hours) */
/* Outputs :
   omega : solar hour angle (in radians) */
/* The procedure "solar_hour_angle_h" supplies an average value of the solar 
   hour angle (in radians) for a whole solar hour, taking into account only the
   portion of the solar hour with the sun standing above the horizon. Returns 0
   if OK, 1 otherwise. */
int solar_hour_angle_h (double phi_g, double delta, double t,
			 double *omega);

/*********************************/
/* SUNRISE, SUNSET AND DAYLENGTH */
/*********************************/

/* Source : */
/* Inputs :
   phi_g       : latitude of site (in radians, positive to North)
   delta       : solar declination angle (in radians)
   gamma_riset : solar elevation near sunrise/sunset:
   - set to  0.0 for astronomical sunrise/sunset
   - set to -1.0 for refraction corrected sunrise/sunset. */
/* Outputs :
   omega_sr : sunrise solar hour angle (in radians)
   omega_ss : sunset solar hour angle (in radians) */
/* The procedure "sunrise_hour_angle" supplies the sunrise and sunset hour 
   angles (in radians). Due to the dimension of the solar disk and the effect 
   of the atmospheric refraction, the edge of the solar disk will just appear 
   (disappear) at the horizon at sunrise (at sunset) when the calculated 
   astronomical elevation is 50'. Returns 0 if OK, 1 otherwise. */
int sunrise_hour_angle (double phi_g, double delta, double gamma_riset,
			 double *omega_sr, double *omega_ss);

/* Same algorithm for sunset only */
double sunset(double phi, double delta);

/* Source : */
/* Inputs :
   omega_sr : sunrise hour angle (in radians)
   omega_ss : sunset hour angle (in radians) */
/* Outputs :
   t_sr : time of astronomical sunrise (in decimal hours)
   t_ss : time of astronomical sunset (in decimal hours)
   S0   : astronomical daylength (in decimal hours) */
/* The procedure "timerise_daylength" supplies the times of astronomical 
   sunrise and sunset, and the astronomical daylength, all in LAT decimal 
   hours. Returns 0 if OK, 1 otherwise. */
int timerise_daylength (double omega_sr, double omega_ss,
			 double *t_sr, double *t_ss, double *S0);

/****************************/
/* CHANGING THE TIME SYSTEM */
/****************************/

/* Source : Gruter (ed.) (1984) */
/* Inputs :
   day_angle   : day angle (in radians)
   lambda      : longitude of the site (in radians, positive to East)
   lambda_ref  : reference longitude of the time zone (in radians)
   summer_corr : correction for summer time (integer hours) */
/* Outputs :
   dt : Offset between local mean time (LMT) and local apparent time (LAT) (in 
   decimal hours) */
/* The procedure "LMT_to_LAT computes the difference (in decimal hours) between
   the LAT (local apparent time) and the LMT (local mean time or clock time) 
   systems at solar noon. Two stages:
   - the first stage calculates the equation of time, ET, wich allows for 
   perturbations in the rotational and angular orbital speed of the Earth.
   - the second stage handles the difference between the longitude of the site 
   under consideration and the reference time zone longitude for the site. A 
   summer time correction must be added for some countries. */
double LMT_to_LAT(double day_angle, double lambda, double lambda_ref, int summer_corr);

/* Source : */
/* Inputs :
   UT          : Universal Time (in decimal hours)
   day_angle   : day angle (in radians)
   lambda      : longitude of the site (in radians, positive to East) */
/* Outputs :
   LAT : local apparent time or solar time or true solar time (TST) (in decimal
   hours) */
/* The procedure "UT_to_LAT computes the conversion of the UT (Universal time) 
   into the LAT (local apparent time) systems at solar noon (in decimal hours).
   First, the equation of time, ET, is computed (in decimal hours), wich allows
   for perturbations in the rotational and angular orbital speed of the Earth.
   Returns 0 if OK, 1 otherwise. */
int UT_to_LAT (double UT, double day_angle, double lambda,
		double *LAT);

/***************************************/
/* POSITION OF THE SUN IN THE SKY FAST */
/***************************************/

// Structure representing the sun position w.r.t. a location on Earth and a day
typedef struct {
    
    angle_t phi;   /* phi    : geocentric latitude of the site, positive to North */
    angle_t delta; /* delta  : solar declination angle */
    
    double sin_phi_sin_delta;
    double cos_phi_cos_delta;
    
} S_GEO_LOCATION_FAST;

void init_geo_location_fast(S_GEO_LOCATION_FAST *p_loc,
                                          const angle_t* phi, const angle_t* delta);
void elevation_sun_fast(const S_GEO_LOCATION_FAST *p_loc, const angle_t* omega,
                                      angle_t *p_gamma);
void elevation_zenith_sun_fast(const S_GEO_LOCATION_FAST *p_loc, const angle_t* omega,
                                             angle_t *p_gamma, angle_t *p_theta);
void azimuth_sun_fast(const S_GEO_LOCATION_FAST *p_loc, const angle_t* omega, const angle_t* gamma,
                                    angle_t *p_alpha);

// Structure representing a tilted plane
typedef struct {
    
    S_GEO_LOCATION_FAST loc; // sun position w.r.t. a location on Earth and a day
    
    angle_t alpha; // azimuth
    angle_t beta;  // tilt
    
    double A;
    double B;
    double C;   

} S_TILTED_PLANE_FAST;

// Create a tilted plane at location p_loc with alpha azimuth and beta tilt 
void init_tilted_plane_fast(S_TILTED_PLANE_FAST *p_tp,
                                          const S_GEO_LOCATION_FAST *p_loc, const angle_t* alpha, const angle_t* beta);
void cos_incident_angle_fast(const S_TILTED_PLANE_FAST *p_tp, const angle_t* omega,
                                           double *p_costhetai);

/**********************************/
/* POSITION OF THE SUN IN THE SKY */
/**********************************/

/* Source : */
/* Inputs :
   phi_g : latitude of site (in radians, positive to North)
   delta : solar declination angle (in radians)
   omega : solar hour angle (in radians) */
/* Outputs :
   gamma : solar altitude angle (in radians)
   theta : solar zenithal angle (in radians) */
/* The procedure "elevation_zenith_sun" computes the solar elevation (or 
   altitude) angle and the solar zenithal (or incidence) angle. These two 
   angles are complementary. Returns 0 if OK, 1 otherwise. */
int elevation_zenith_sun (double phi_g, double delta, double omega,
			   double *gamma, double *theta);

/* Source : */
/* Inputs :
   phi_g : latitude of site (in radians, positive to North)
   delta : solar declination angle (in radians)
   omega : solar hour angle (in radians)
   gamma : solar altitude angle (in radians) */
/* Outputs :
   alpha : solar azimuthal angle (in radians) */
/* The procedure "azimuth_sun" computes the solar azimuth angle in the Northern
   hemisphere. The azimuth angle has a positive value when the sun is to the 
   west of South, i.e. during the afternoon in solar time. For the Southern 
   hemisphere, the azimuth angle is measured from North. Returns 0 if OK, 1 
   otherwise. */
int azimuth_sun (double phi_g, double delta, double omega, double gamma,
		  double *alpha);

/* Source : Lucien Wald */
/* Inputs :
   phi      : latitude of site (in radians, positive to North)
   alpha    : solar azimuthal angle using Solar Geometry convention (in radians) */
/* Outputs :
   alpha_iso: solar azimuthal angle using ISO convention (in radians) */
/* The procedure "azimuth_sg_to_azimuth_iso" converts the solar azimuthal angle expressed
   using Solar Geometry convention (in the Northern hemisphere, the azimuth angle 
   has a positive value when the sun is to the west of South, i.e. during the afternoon 
   in solar time ; for the Southern hemisphere, the azimuth angle is measured from North) 
   into the solar azimuthal angle expressed using ISO convention (0 = North, East is positive). */
double azimuth_sg_to_azimuth_iso(double phi, double alpha);

/* Source : Lucien Wald */
/* Inputs :
   phi      : latitude of site (in radians, positive to North)
   alpha_iso: solar azimuthal angle using ISO convention (in radians) */
/* Outputs :
   alpha    : solar azimuthal angle using Solar Geometry convention (in radians) */
/* The procedure "azimuth_iso_to_azimuth_sg" converts the solar azimuthal angle expressed
   using ISO convention (0 = North, East is positive) into the solar azimuthal angle expressed
   using Solar Geometry convention (in the Northern hemisphere, the azimuth angle 
   has a positive value when the sun is to the west of South, i.e. during the afternoon 
   in solar time ; for the Southern hemisphere, the azimuth angle is measured from North). */
double azimuth_iso_to_azimuth_sg(double phi, double alpha_iso);

/* Source : */
/* Inputs :
   phi_g : latitude of site (in radians, positive to North)
   delta : solar declination angle (in radians)
   omega : solar hour angle (in radians)
   gamma : solar altitude angle (in radians)
   beta  : plane tilt angle (in radians)
   alpha : plane azimuthal angle (in radians) */
/* Outputs :
   costhetai : cosinus of the incidence angle of the sun wrt plane */
/* The procedure "cos_incident_angle" computes the cosinus of the incidence angle 
   of the sun wrt a tilted plane. Returns 0 if OK, 1 otherwise. */
int cos_incident_angle(double phi_g, double delta, double omega, double gamma,
        double beta, double alpha,
        double *costhetai);

/********************************/
/* EXTRATERRESTRIAL IRRADIATION */
/********************************/

/* Source : Gruter (ed.) (1984) */
/* Inputs :
   day_angle : day angle (in radians) */
/* Outputs :
   eccentricity : correction for Earth orbit eccentricity */
/* The procedure "corr_distance" computes the correction for the variation of 
   sun-earth distance from its mean value (also known as eccentricity). It is a
   fucntion of time, but a single (average) value per day is enough for 
   practical calculations. */
double corr_distance(double day_angle);

/* Source : */
/* Inputs :
   IOj : extraterrestrial solar irradiance normal to beam for day j; I0j=I0*fj
   theta : solar incidence angle or solar zenithal angle */
/* Outputs :
   G0 : extraterrestrial global solar irradiation (in Wh/m2) */
/* The procedure "G0_normal" delivers the extraterrestrial solar irradiance 
   normal to beam for day j. Returns 0 if OK, 1 otherwise. */
int G0_normal (double I0j, double theta, double *G0);

/* Source : */
/* Inputs :
   phi_g        : latitude of site (in radians, positive to North)
   eccentricity : correction for Earth orbit eccentricity
   delta        : solar declination angle (in radians)
   omega1       : solar hour angle at beginning of the period (in radians)
   omega2       : solar hour angle at end of the period (in radians) */
/* Outputs :
   G0_12 : extraterrestrial solar irradiation (in Wh/m2) */
/* The procedure "G0_general" delivers the extraterrestrial solar irradiation 
   incident on an horizontal surface in the general case (in Wh/m2). Returns 0
   if OK, 1 otherwise */
int G0_general (double phi_g, double eccentricity, double delta, double omega1, double omega2,
		 double *G0_12);

/* Source : */
/* Inputs :
   phi_g        : latitude of site (in radians, positive to North)
   eccentricity : correction for Earth orbit eccentricity
   delta        : solar declination angle (in radians) */
/* Outputs :
   G0d : daily extraterrestrial solar irradiation (in Wh/m2) */
/* The procedure "G0_day" delivers the extraterrestrial solar irradiation 
   incident on an horizontal surface in case of daily values (in Wh/m2), i.e.
   omega1 = omega_sr = -omega_ss  et omega2 = omega_ss. Returns 0 if OK, 1 
   otherwise.
   REMARK: It is a special case of G0_general with the sunrise and sunset 
   angles as integration limits. */
int G0_day (double phi_g, double eccentricity, double delta,
	     double *G0d);

/* Source : */
/* Inputs :
   phi_g        : latitude of site (in radians, positive to North)
   eccentricity : correction for Earth orbit eccentricity
   delta        : solar declination (in radians) */
/* Outputs :
   G0h[1..24] : 24 hourly extraterrestrial solar irradiation (in Wh/m2) */
/* The procedure "G0_hours_profile" delivers the extraterrestrial solar 
   irradiation incident on an horizontal surface in case of hourly values, for
   the 24 integral hours in a given day (in Wh/m2), i.e. |omega1-omega2| = 
   Pi/12. Returns 0 if OK, 1 otherwise */
int G0_hours_profile(double phi_g, double eccentricity, double delta,
        double *G0h);

/* Source : */
/* Inputs :
   phi_g        : latitude of site (in radians, positive to North)
   eccentricity : correction for Earth orbit eccentricity
   delta        : solar declination (in radians)
   t            : solar time i.e. LAT (0..24 decimal hours) */
/* Outputs :
   G0h : hourly extraterrestrial solar irradiation (in Wh/m2) */
/* The procedure "G0_hour" delivers the extraterrestrial solar irradiation 
   incident on an horizontal surface for a specific hour in a given day (in 
   Wh/m2), i.e. |omega1-omega2| = Pi/12. t is taken as the mid hour for 
   computation of the hourly value of extraterrestrial solar irradiation. 
   Returns 0 if OK, 1 otherwise */
int G0_hour (double phi_g, double eccentricity, double delta, double t,
	      double *G0h);

/***********************************************/
/* MONTHLY AVERAGES OF SOLAR INPUTS PARAMETERS */
/***********************************************/
/* Source : */
/* Inputs :
   month_number : month number (1..12)
   year_number  : year number (4 digits)
   phi_g        : latitude of site (in radians, positive to North)
   lambda       : longitude of site (in radians, positive to East)
   gamma_riset  : solar elevation near sunrise/sunset:
   - set to  0.0 for astronomical sunrise/sunset
   - set to -1.0 for refraction corrected sunrise/sunset. */
/* Outputs :        monthly average of...
   day_angle_m    : ... day angle (in radians)
   delta_m        : ... solar declination angle (in radians)
   omega_ss_m     : ... sunset hour angle (in radians)
   S0_m           : ... astronomical daylength (in decimal hours)
   eccentricity_m : ... eccentricity
   G0d_m          : ... daily extraterrestrial irradiation (in Wh/m2)
   G0h_m[1..24]   : ... 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
/* The procedure "monthly_averages" computes directly the monthly average 
   values of solar parameters : day angle (in radians), eccentricity, 
   declination (in radians), sunset hour angle (in radians), daylength (in 
   decimal hours), daily extraterrestrial irradiation (in Wh/m2) and the 24 
   hourly extraterrestrial solar irradiation (in Wh/m2) . Returns 0 if OK, 1 
   otherwise */
int monthly_averages (int month_number, int year_number,
               double phi_g, double lambda, double gamma_riset, 
               double *day_angle_m, double *delta_m, double *omega_ss_m,
		       double *S0_m, double *eccentricity_m,
		       double *G0d_m, double *G0h_m);

/****************************************************/
/* YEARLY AVERAGES OF MONTHLY SOLAR PARAMETERS      */
/* (LONG TERM MEANS OF MONTHLY MEANS OF DAILY SUMS) */
/****************************************************/
/* Source : */
/* Inputs :
   month_number : month number (1..12)
   year_start   : starting year of the considered period (4 digits)
   year_end     : ending year of the considered period (4 digits)
   phi_g        : latitude of site (in radians, positive to North)
   lambda       : longitude of site (in radians, positive to East)
   gamma_riset  : solar elevation near sunrise/sunset:
   - set to  0.0 for astronomical sunrise/sunset
   - set to -1.0 for refraction corrected sunrise/sunset. */
/* Outputs :        yearly average of...
   day_angle_y    : ... day angle (in radians)
   delta_y        : ... solar declination angle (in radians)
   omega_ss_y     : ... sunset hour angle (in radians)
   S0_y           : ... astronomical daylength (in decimal hours)
   eccentricity_y : ... eccentricity
   G0d_y          : ... daily extraterrestrial irradiation (in Wh/m2)
   G0h_y[1..24]   : ... 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
/* The procedure "yearly_averages" computes directly the yearly average over a
   defined period of years of monthly average values of solar parameters : day 
   angle (in radians), eccentricity, declination (in radians), sunset hour 
   angle (in radians), daylength (in decimal hours), daily extraterrestrial 
   irradiation (in Wh/m2) and 24 hourly extraterrestrial solar irradiation 
   (in Wh/m2). Returns 0 if OK, 1 otherwise */
int yearly_averages(int month_number, int year_start, int year_end,
              double phi_g, double lambda, double gamma_riset, 
              double *day_angle_y, double *delta_y, double *omega_ss_y, 
		      double *S0_y, double *eccentricity_y,
		      double *G0d_y, double *G0h_y);

/*********************************************/
/* SOLAR INPUTS PARAMETERS FOR A CERTAIN DAY */
/*********************************************/
/* Source : */
/* Inputs :
   day_of_month : day of the month (1..31)
   month_number : month number (1..12)
   year_number  : year number (4 digits)
   phi_g        : latitude of site (in radians, positive to North)
   lambda       : longitude of site (in radians, positive to East)
   gamma_riset  : solar elevation near sunrise/sunset:
   - set to  0.0 for astronomical sunrise/sunset
   - set to -1.0 for refraction corrected sunrise/sunset. */
/* Outputs :
   day_angle    : day angle (in radians)
   delta        : solar declination angle (in radians)
   omega_ss     : sunset hour angle (in radians)
   S0           : astronomical daylength (in decimal hours)
   eccentricity : eccentricity
   G0d          : daily extraterrestrial irradiation (in Wh/m2)
   G0h[1..24]   : 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
/* The procedure "solar_parameters_day" computes the solar geometry related 
   values for a certain day : day angle (in radians), eccentricity, 
   declination (in radians), sunset hour angle (in radians), daylength (in 
   decimal hours), daily extraterrestrial irradiation (in Wh/m2) and the 24 
   hourly extraterrestrial solar irradiation (in Wh/m2). Returns 0 if OK, 1 
   otherwise. 
   REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar. */
int solar_parameters_day(int day_of_month, int month_number, int year_number,
                double phi_g, double lambda, double gamma_riset, 
                double *day_angle, double *delta, double *omega_ss, 
                double *S0, double *eccentricity, 
                double *G0d, double *G0h);

/******************************************************************/
/* SOLAR INPUTS PARAMETERS FOR MONTHLY AVERAGE IRRADIATION MODELS */
/******************************************************************/
/* Source : */
/* Inputs :
   month_number : month number (1..12)
   phi_g        : latitude of site (in radians, positive to North)
   gamma_riset  : solar elevation near sunrise/sunset:
   - set to  0.0 for astronomical sunrise/sunset
   - set to -1.0 for refraction corrected sunrise/sunset. */
/* Outputs :        average ... for the given month
   day_angle_avg    : day angle (in radians)
   delta_avg        : solar declination angle (in radians)
   omega_ss_avg     : sunset hour angle (in radians)
   S0_avg           : astronomical daylength (in decimal hours)
   eccentricity_avg : eccentricity
   G0d_avg          : daily extraterrestrial irradiation (in Wh/m2)
   G0h_avg[1..24]   : 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
/* The procedure "solar_parameters_acg" computes the solar geometry related 
   values for monthly average irradiation models : day angle (in radians), 
   eccentricity, declination (in radians), sunset hour angle (in radians), 
   daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2) 
   and the 24 hourly extraterrestrial solar irradiation (in Wh/m2). Returns 0 
   if OK, 1 otherwise. 
   REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar. */
int solar_parameters_avg (int month_number,
			   double phi_g, double gamma_riset,
			   double *day_angle_avg,
			   double *delta_avg,
			   double *omega_ss_avg, double *S0_avg,
			   double *eccentricity_avg,
			   double *G0d_avg, double *G0h_avg);

/******************************************************************/
/* SOLAR INPUTS PARAMETERS FOR MONTHLY MAXIMUM IRRADIATION MODELS */
/******************************************************************/
/* Source : */
/* Inputs :
   month_number : month number (1..12)
   phi_g        : latitude of site (in radians, positive to North)
   gamma_riset  : solar elevation near sunrise/sunset:
   - set to  0.0 for astronomical sunrise/sunset
   - set to -1.0 for refraction corrected sunrise/sunset. */
/* Outputs :        average ... for the given month
   day_angle_max    : day angle (in radians)
   delta_max        : solar declination angle (in radians)
   omega_ss_max     : sunset hour angle (in radians)
   S0_max           : astronomical daylength (in decimal hours)
   eccentricity_max : eccentricity
   G0d_max          : daily extraterrestrial irradiation (in Wh/m2)
   G0h_max[1..24]   : 24 hourly extraterrestrial solar irradiation (in Wh/m2)
 */
/* The procedure "solar_parameters_acg" computes the solar geometry related 
   values for monthly average irradiation models : day angle (in radians), 
   eccentricity, declination (in radians), sunset hour angle (in radians), 
   daylength (in decimal hours), daily extraterrestrial irradiation (in Wh/m2) 
   and the 24 hourly extraterrestrial solar irradiation (in Wh/m2). Returns 0 
   if OK, 1 otherwise. 
   REMARK: gamma_riset set to 0.0 in the original procedure by Aguiar. */
int solar_parameters_max (int month_number,
			   double phi_g, double gamma_riset,
			   double *day_angle_max,
			   double *delta_max,
			   double *omega_ss_max, double *S0_max,
			   double *eccentricity_max,
			   double *G0d_max, double *G0h_max);

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
int intervals_omega_tilted_plane (double phi_g, double delta, double omega_ss, double beta, double alpha,
                double *v_om, int *p_nb);

#ifdef	__cplusplus
}
#endif

#endif /* __H_solar_geometry */
