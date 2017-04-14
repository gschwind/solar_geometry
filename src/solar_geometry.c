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

#define __C_solar_geometry

#include "solar_geometry.h"

#include <math.h>
#include <stdio.h>
#include <assert.h>


/*************/
/* CONSTANTS */
/*************/

// #days of each month (0-based) for NON leap year
static const int NB_DAYS_OF_MONTH[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
// #days of each month (0-based) for leap year
static const int NB_DAYS_OF_MONTH_BISSEXTILE_YEAR[12] = {31,29,31,30,31,30,31,31,30,31,30,31};

// 1-based index of the last day of each month (1-based) for NON leap year
static const int FIRST_DAY_OF_MONTH[13]                 = {0,31,59,90,120,151,181,212,243,273,304,334,365};
// 1-based index of the last day of each month (1-based) for leap year
static const int FIRST_DAY_OF_MONTH_BISSEXTILE_YEAR[13] = {0,31,60,91,121,152,182,213,244,274,305,335,366};

// name of each month (1-based)
static const char NAME_OF_MONTH[][4] = { "", "jan", "feb", "mar", "apr", "may", "jun", 
                                       "jul", "aug", "sep", "oct", "nov", "dec" };


#define bissextile(yyyy) ((((yyyy)%4)==0&&((yyyy)%100)!=0)||((yyyy)%400)==0)

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))

#ifndef M_PI
#define M_PI 3.141592653589793
#define M_PI_2 1.570796326794897
#endif
#define M_PI_3 1.047197551196598
#define M_PI_12 0.261799387799149
#define M_PI_24 0.130899693899575

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

int make_julian_day(int day_of_month, int month_number, int year_number,
  int *julian_day)
{
  int ier, julien;

  assert((day_of_month > 0) && (day_of_month < 32) && (month_number > 0)
      && (month_number < 13) && (year_number > 0));

  ier = 1;
  if ((day_of_month > 0) && (day_of_month < 32) && (month_number > 0)
      && (month_number < 13) && (year_number > 0))
    {
      ier = 0;
      julien = day_of_month + FIRST_DAY_OF_MONTH[month_number - 1];
      if (bissextile(year_number) && (month_number > 2))
        /*
         * February of leap year
         */
        julien = julien + 1;
      *julian_day = julien;
    }
  return (ier);
}

int julian_to_date(int year_number, int julian_day,
    int *day_of_month, int *month_number)
{
  const int *tab;
  int ier, m, jmax = 365;

  assert((julian_day > 0) && (julian_day <= jmax) && (year_number > 0));

  ier = 1;
  if (bissextile(year_number))
    { /* leap year */
      jmax = jmax + 1;
      tab = FIRST_DAY_OF_MONTH_BISSEXTILE_YEAR;
    }
  else
    {
      tab = FIRST_DAY_OF_MONTH;
    }

  if ((julian_day > 0) && (julian_day <= jmax) && (year_number > 0))
    {
      ier = 0;
      for (m = 0; m < 12; m++)
        {
          if ((julian_day > tab[m]) && (julian_day <= tab[m + 1]))
            {
              *month_number = m + 1;
              *day_of_month = julian_day - tab[m];
              break;
            }
          else if (julian_day > tab[11])
            {
              *month_number = 12;
              *day_of_month = julian_day - tab[11];
              break;
            }
        }
    }

  return (ier);
}

inline int nbdays_month(int year_number, int month_number)
{
    assert((year_number > 0) && (month_number > 0) && (month_number < 13));
    
	if (bissextile(year_number)) {
		return NB_DAYS_OF_MONTH_BISSEXTILE_YEAR[month_number - 1];
	} else {
		return NB_DAYS_OF_MONTH[month_number - 1];
	}
}

int number_to_name_month(int month_number, char *month_name)
{
  int ier;

  assert((month_number > 0) && (month_number < 13));
  
  ier = 1;
  if ((month_number > 0) && (month_number < 13))
    {
      ier = 0;
      sprintf(month_name, "%s", NAME_OF_MONTH[month_number]);
    }

  return (ier);
}

inline extern double Day_Angle(int julian_day)
{
  assert((julian_day > 0) && (julian_day <= 366));

  return (double)julian_day * 2.0 * Pi / 365.2422;
}

inline extern double declination_sun(int year_number, int julian_day, double lambda)
{
	double n0, t1, wt;

	double const b1 = 0.0064979;
	double const b2 = 0.4059059;
	double const b3 = 0.0020054;
	double const b4 = -0.0029880;
	double const b5 = -0.0132296;
	double const b6 = 0.0063809;
	double const b7 = 0.0003508;
	double const w0 = 2.0 * Pi / 365.2422;

    assert ((julian_day > 0) && (julian_day <= 366));

    /*
    * n0 : spring-equinox time expressed in days from the beginning of the year i.e.
    * the time in decimal days elapsing from 00:00 hours Jan 1st to the spring equinox
    * at Greenwich in a given year.
	* ((year_number - 1957) >> 2)  <=> INT[(year_number - 1957)/4] where INT is the integral value.
    * t1 : time in days, from the spring equinox.
    * 0.5 represents the decimal day number at noon on Jan 1st at Greenwich.
    */
	n0 = 78.8946 + 0.2422 * (year_number - 1957) - ((year_number - 1957) >> 2);
    t1 = -0.5 - lambda / (2 * Pi) - n0;
	wt = w0 * (julian_day + t1);

	return b1 + b2 * sin(wt) + b3 * sin(2 * wt) + b4 * sin(3 * wt)
			  + b5 * cos(wt) + b6 * cos(2 * wt) + b7 * cos(3 * wt);
}

int declination_sun_month(int month_number, int type_use, 
  double *delta_month)
{
  const double deg_rad = (Pi / 180.0); /* converts decimal degrees into radians */
  int tab_julian_day[12] =
    { 17, 46, 75, 105, 135, 162, 198, 228, 259, 289, 319, 345 };
  int tab_julian_day_max[12] =
    { 29, 57, 89, 119, 150, 173, 186, 217, 248, 278, 309, 339 };
  int ier, julian_day = 0;
  double day_angle, jm, c1, c2, c3, c4;

  assert((type_use >= 0) && (type_use < 2));

  ier = 1;
  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 0)
        julian_day = tab_julian_day[month_number - 1];
      if (type_use == 1)
        julian_day = tab_julian_day_max[month_number - 1];
      ier = 0;
    }

  day_angle = Day_Angle(julian_day);
  
  if (ier != 0)
    return (ier);

  jm = day_angle;
  c1 = 0.3978;
  c2 = 80.2 * deg_rad; /* 1.4000 in SSA manual */
  c3 = 1.92 * deg_rad; /* 0.0355 in SSA manual */
  c4 = 2.80 * deg_rad; /* 0.0489 in SSA manual */

  *delta_month = asin(c1 * sin(jm - c2 + c3 * sin(jm - c4)));
  return (ier);
}

inline extern double solar_hour_angle (double t)
{
//  assert((t >= 0.0) && (t <= 24.0));

  return (t - 12.0) * Pi / 12.0;
}

inline extern double omega_to_LAT (double omega)
{
  assert((omega >= -Pi) && (omega <= Pi));

  return 12.0 * (1.0 + omega / Pi);
}

inline extern double geogr_to_geoce(double phi_g)
{
    double const CC = 0.99330552; /* Correction factor for converting geographic latitude
                                   * into geocentric latitude. 
                                   * CC=(Rpole/Requator)**2
                                   * Rpole=6356.752, Requator=6378.137
                                   */
                                   
    assert((phi_g >= -(Pi / 2.0 - 0.0002)) || (phi_g <= (Pi / 2.0 - 0.0002)));
	
    return atan(tan(phi_g) * CC);
}

int solar_hour_angle_h(double phi_g, double delta, double t,
  double *omega)
{
  int ier;
  double omega_sr, omega_ss, omega1, omega2;

  ier = 1;
  ier = sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);

  omega1 = (t - 1.0 - 12.0) * Pi / 12.0;
  if (omega1 < omega_sr)
    omega1 = omega_sr;

  omega2 = (t - 12.0) * Pi / 12.0;
  if (omega2 > omega_ss)
    omega2 = omega_ss;

  *omega = (omega1 + omega2) / 2.0;

  return (ier);
}

/*********************************/
/* SUNRISE, SUNSET AND DAYLENGTH */
/*********************************/

int sunrise_hour_angle(double phi_g, double delta, double gamma_riset,
    double *omega_sr, double *omega_ss)
{
    static double deg_rad = (Pi / 180.0); /* converts decimal degrees into radians */
    int ier;
    double horizon, max_delta, cos_omega_sunset = 1, omega_sunset = 0;
    double phi;

    assert((gamma_riset == 0.0) || (gamma_riset == -1.0));

    ier = 1;
    if ((gamma_riset == 0.0) || (gamma_riset == -1.0))
        ier = 0;
    horizon = (-50.0 / 60.0) * deg_rad; /* horizon, -50' in radians */
    if (gamma_riset >= horizon)
        horizon = gamma_riset;

    phi = geogr_to_geoce(phi_g);
    max_delta = 23.45 * deg_rad;

    assert((fabs(phi) < (Pi / 2.0)) && (fabs(delta) <= max_delta));

    if ((fabs(phi) < (Pi / 2.0)) && (fabs(delta) <= max_delta) && (ier == 0))
    {
        cos_omega_sunset = (sin(horizon) - (sin(phi) * sin(delta))) / (cos(phi) * cos(delta));
        ier = 0;
    }
    else
        ier = 1;

    if (cos_omega_sunset >= 1.0) /* the sun is always below the horizon: polar night */
        omega_sunset = 0.0;
    else if (cos_omega_sunset <= -1.0) /* the sun is always above the horizon: polar day */
        omega_sunset = Pi;
    else
        omega_sunset = acos(cos_omega_sunset);

    *omega_sr = -omega_sunset;
    *omega_ss = omega_sunset;
    return (ier);
}

inline double sunset(double phi, double delta)
{
    //double cos_omega_sunset = (sin(0.0) - (sin(phi) * sin(delta))) / (cos(phi) * cos(delta));
	double cos_omega_sunset = - tan(phi) * tan(delta);
    if (cos_omega_sunset >= 1.0) /* the sun is always below the horizon: polar night */
		return 0.0;
    else if (cos_omega_sunset <= -1.0) /* the sun is always above the horizon: polar day */
		return Pi;
    else
	    return acos(cos_omega_sunset);
}

int timerise_daylength(double omega_sr, double omega_ss, 
  double *t_sr, double *t_ss, double *S0)
{
  int ier;

  assert((omega_sr >= -Pi) && (omega_sr <= 0.0) && (omega_ss >= 0.0) && (omega_ss <= Pi));
  
  ier = 1;
  if ((omega_sr >= -Pi) && (omega_sr <= 0.0) && (omega_ss >= 0.0) && (omega_ss <= Pi))
    {
      ier = 0;
      /*
       * alternative way
       */
      /*
       * ier = omega_to_LAT(omega_sr,&t_sr); 
       * if(ier == 0) ier = omega_to_LAT(omega_ss,&t_ss); 
       * if(ier != 0) return(ier);
       */
      *t_sr = 12.0 + omega_sr * 12.0 / Pi;
      *t_ss = 12.0 + omega_ss * 12.0 / Pi;
      *S0 = *t_ss - *t_sr;
    }
  return (ier);
}

/****************************/
/* CHANGING THE TIME SYSTEM */
/****************************/

inline double LMT_to_LAT(double day_angle, double lambda, double lambda_ref, int summer_corr)
{
    const double deg_rad = (Pi / 180.0); /* converts decimal degrees into radians */
    double const a1 = -0.128;
    double const a2 = -0.165;
    double const a3 = 2.80 * deg_rad;
    double const a4 = 19.70 * deg_rad;
    double ET;

    assert((day_angle > 0.0) && (day_angle < (2.0 * Pi * 1.0021)) 
        && (fabs(lambda) <= Pi) && (fabs(lambda_ref) <= Pi));

    ET = a1 * sin(day_angle - a3) + a2 * sin(2.0 * day_angle + a4);
    return ET + ((lambda - lambda_ref) * 12.0 / Pi) - (double) summer_corr;
}

int UT_to_LAT(double UT, double day_angle, double lambda, 
  double *LAT)
{
  const double deg_rad = (Pi / 180.0); /* converts decimal degrees into radians */
  int ier;
  double a1, a2, a3, a4, ET;

  assert((UT >= 0.0) && (UT <= 24.0) 
      && (day_angle > 0.0) && (day_angle < (2.0 * Pi * 1.0021)) 
      && (fabs(lambda) <= Pi));
  
  ier = 1;
  a1 = -0.128;
  a2 = -0.165;
  a3 = 2.80 * deg_rad;
  a4 = 19.70 * deg_rad;
  if ((UT >= 0.0) && (UT <= 24.0) 
   && (day_angle > 0.0) && (day_angle < (2.0 * Pi * 1.0021)) 
   && (fabs(lambda) <= Pi))
    {
      ier = 0;
      ET = a1 * sin(day_angle - a3) + a2 * sin(2.0 * day_angle + a4);
      *LAT = UT + ET + (lambda * 12.0 / Pi);
      if (*LAT < 0)
        *LAT += 24.0;
      if (*LAT > 24.0)
        *LAT -= 24.0;
    }

  return (ier);
}

/***************************************/
/* POSITION OF THE SUN IN THE SKY FAST */
/***************************************/

void init_geo_location_fast(S_GEO_LOCATION_FAST *p_loc, 
                            const angle_t* phi, const angle_t* delta)
{	
	double sin_phi;
	double cos_phi;
	double sin_delta;
	double cos_delta;

	assert(p_loc != NULL);
    assert(phi != NULL);
    assert(delta != NULL);

	sin_phi   = get_sin(phi);
	cos_phi   = get_cos(phi);
	sin_delta = get_sin(delta);
	cos_delta = get_cos(delta);

	p_loc->sin_phi_sin_delta = sin_phi * sin_delta;
	p_loc->cos_phi_cos_delta = cos_phi * cos_delta;
    
    // copy angle structures once sinus and cosinus are computed
    p_loc->phi   = *phi;
    p_loc->delta = *delta;
}

void elevation_sun_fast(const S_GEO_LOCATION_FAST *p_loc, const angle_t* omega, 
                        angle_t *p_gamma)
{

	double cos_omega;
    double sin_gamma;

	assert(p_loc != NULL);
    assert(omega != NULL);
    assert(p_gamma != NULL);
    
	cos_omega = get_cos(omega);
	sin_gamma = p_loc->sin_phi_sin_delta + p_loc->cos_phi_cos_delta*cos_omega;

    if (sin_gamma < 0.0)
        sin_gamma = 0.0; // force gamma sun = 0 before sunrise or after sunset
    *p_gamma = init_angle_sin(asin(sin_gamma), sin_gamma);

}

void elevation_zenith_sun_fast(const S_GEO_LOCATION_FAST *p_loc, const angle_t* omega, 
                               angle_t *p_gamma, angle_t *p_theta)
{
	double cos_omega;
    double sin_gamma;

    assert(p_loc != NULL);
    assert(omega != NULL);
    assert(p_gamma != NULL);
    assert(p_theta != NULL);
    
	cos_omega = get_cos(omega);

    sin_gamma = p_loc->sin_phi_sin_delta + p_loc->cos_phi_cos_delta*cos_omega;
    if (sin_gamma < 0.0)
        sin_gamma = 0.0; // force gamma sun = 0 before sunrise or after sunset
    *p_gamma = init_angle_sin(asin(sin_gamma), sin_gamma);
    
	*p_theta = pi_2_minus_angle(p_gamma); // PI/2 - gamma

	
}

void azimuth_sun_fast(const S_GEO_LOCATION_FAST *p_loc, const angle_t* omega, const angle_t* gamma, 
                      angle_t *p_alpha)
{

	double sin_gamma;
	double cos_gamma;
	double phi;
	double sin_phi;
	double cos_phi;
	double sin_delta;
	double cos_delta;
	double sin_omega;
	double cos_as; // azimuth sun
	double sin_as;
	double as;

    assert(p_loc != NULL);
    assert(omega != NULL);
    assert(gamma != NULL);
    assert(p_alpha != NULL);
    
	sin_gamma = get_sin(gamma);
	cos_gamma = get_cos(gamma);
	phi       = get_angle(&p_loc->phi);
	sin_phi   = get_sin(&p_loc->phi);
	cos_phi   = get_cos(&p_loc->phi);
	sin_delta = get_sin(&p_loc->delta);
	cos_delta = get_cos(&p_loc->delta);
	sin_omega = get_sin(omega);

	cos_as = (sin_phi * sin_gamma - sin_delta) / (cos_phi * cos_gamma); // azimuth sun
	if (phi < 0.0)
		cos_as = -cos_as; /* Southern hemisphere */
		
	sin_as = cos_delta * sin_omega / cos_gamma;

	if (cos_as > 1.0)
		cos_as = 1.0;
	else if (cos_as < -1.0)
		cos_as = -1.0;

	as = acos(cos_as);
	if (sin_as >= 0.0)
		*p_alpha = init_angle_cos_sin( as, cos_as, sin_as);
	else
		*p_alpha = init_angle_cos_sin(-as, cos_as, sin_as);

}

void init_tilted_plane_fast(S_TILTED_PLANE_FAST *p_tp, 
                            const S_GEO_LOCATION_FAST *p_loc, const angle_t* alpha, const angle_t* beta)
{

	double sin_alpha;	
	double cos_alpha;	
	double sin_beta;	
	double cos_beta;	
	double phi;
	double sin_phi;
	double cos_phi;
	double sin_delta;
	double cos_delta;

    assert(p_tp != NULL);
    assert(p_loc != NULL);
    assert(alpha != NULL);
    assert(beta != NULL);
    
	sin_alpha = get_sin(alpha);	
	cos_alpha = get_cos(alpha);	
	sin_beta  = get_sin(beta);	
	cos_beta  = get_cos(beta);	
	phi       = get_angle(&p_loc->phi);
	sin_phi   = get_sin(&p_loc->phi);
	cos_phi   = get_cos(&p_loc->phi);
	sin_delta = get_sin(&p_loc->delta);
	cos_delta = get_cos(&p_loc->delta);

	if (phi >= 0.0) {
		p_tp->A = cos_delta * (cos_phi * cos_beta + sin_phi * sin_beta * cos_alpha);
		p_tp->B = cos_delta * sin_beta * sin_alpha;
		p_tp->C = sin_delta * (sin_phi * cos_beta - cos_phi * sin_beta * cos_alpha);
	} else {
		p_tp->A = cos_delta * (cos_phi * cos_beta - sin_phi * sin_beta * cos_alpha);
		p_tp->B = cos_delta * sin_beta * sin_alpha;
		p_tp->C = sin_delta * (sin_phi * cos_beta + cos_phi * sin_beta * cos_alpha);
	}

    // copy angle structures once sinus and cosinus are computed
    p_tp->loc = *p_loc;
    p_tp->alpha = *alpha;
    p_tp->beta  = *beta;
}

void cos_incident_angle_fast(const S_TILTED_PLANE_FAST *p_tp, const angle_t* omega, 
                             double *p_costhetai)
{
	double sin_omega = get_sin(omega);
	double cos_omega = get_cos(omega);

	*p_costhetai = p_tp->A*cos_omega + p_tp->B*sin_omega + p_tp->C;
}

/**********************************/
/* POSITION OF THE SUN IN THE SKY */
/**********************************/

int elevation_zenith_sun(double phi_g, double delta, double omega, 
  double *gamma, double *theta)
{
  int ier;
  double omega_sr, omega_ss;
  double phi;

  ier = 1;
  phi = geogr_to_geoce(phi_g);
  ier = sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  if ((omega < omega_sr) || (omega > omega_ss))
    *gamma = 0.0; // force gamma sun = 0 before sunrise or after sunset
  else
    *gamma = asin(sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega));
  if (*gamma < 0.0)
    *gamma = 0.0;

  *theta = (Pi / 2.0) - *gamma;

  return (ier);
}

int azimuth_sun(double phi_g, double delta, double omega, double gamma,
    double *alpha)
{
  int ier;
  double cos_as, sin_as, as;
  double phi;

  ier = 0;
  phi = geogr_to_geoce(phi_g);
  cos_as = (sin(phi) * sin(gamma) - sin(delta)) / (cos(phi) * cos(gamma));
  if (phi < 0.0)
    cos_as = -cos_as; /* Southern hemisphere */
  sin_as = cos(delta) * sin(omega) / cos(gamma);

  if (cos_as > 1.0)
    cos_as = 1.0;
  else if (cos_as < -1.0)
    cos_as = -1.0;

  as = acos(cos_as);
  if (fabs(as) > Pi)
    ier = 1;
  if (sin_as >= 0.0)
    *alpha = as;
  else
    *alpha = -as;

  return (ier);
}

double azimuth_sg_to_azimuth_iso(double phi, double alpha)
{
    double alpha_iso;
    
    assert(alpha >= -Pi && alpha <= Pi);
    
    // If Northern hemisphere
    if (phi >= 0) {
        alpha_iso = alpha + Pi;
    } else { 
        alpha_iso = -alpha;
        if (alpha_iso < 0.0)
            alpha_iso += 2*Pi;
    }
    
    assert(alpha_iso >= 0.0 && alpha_iso <= 2*Pi);
    return alpha_iso;
}

double azimuth_iso_to_azimuth_sg(double phi, double alpha_iso)
{
    double alpha;
    
    assert(alpha_iso >= 0.0 && alpha_iso <= 2*Pi);

    // If Northern hemisphere
    if (phi >= 0) {
        alpha = alpha_iso - Pi;
    } else {
        alpha = -alpha_iso;
        if (alpha < -Pi)
            alpha += 2*Pi;
    }
    
    assert(alpha >= -Pi && alpha <= Pi);
    return alpha;
}

int cos_incident_angle(double phi_g, double delta, double omega, double gamma, 
    double beta, double alpha,
    double *costhetai)
{
    int ier = 0;
    double sinBeta;
    double cosBeta;
    double A;
    double B;
    double C;
    double phi;
    
    sinBeta = sin(beta);
    cosBeta = cos(beta);

    phi = geogr_to_geoce(phi_g);

    if (phi >= 0.0) {
        A = cos(delta) * (cos(phi) * cosBeta + sin(phi) * sinBeta * cos(alpha));
        B = cos(delta) * sinBeta * sin(alpha);
        C = sin(delta) * (sin(phi) * cosBeta - cos(phi) * sinBeta * cos(alpha));
    } else {
        A = cos(delta) * (cos(phi) * cosBeta - sin(phi) * sinBeta * cos(alpha));
        B = cos(delta) * sinBeta * sin(alpha);
        C = sin(delta) * (sin(phi) * cosBeta + cos(phi) * sinBeta * cos(alpha));
    }

    *costhetai = A * cos(omega) + B * sin(omega) + C;

    return (ier);
}

/********************************/
/* EXTRATERRESTRIAL IRRADIATION */
/********************************/

inline double corr_distance(double day_angle)
{
    const double deg_rad = (Pi / 180.0); /* converts decimal degrees into radians */
    const double a = 2.80 * deg_rad;

    assert((day_angle >= 0.0) && (day_angle <= (2.0 * Pi * 1.0021)));

	return 1.0 + 0.03344 * cos(day_angle - a);
}

int G0_normal(double I0j, double theta, double *G0)
{
  *G0 = I0j * cos(theta);
  return 0;
}

int G0_general(double phi_g, double eccentricity, double delta, double omega1, double omega2, 
    double *G0_12)
{
  int ier;
  double omega_sr, omega_ss, a, b1, b2, c;
  double phi;

  assert(omega2 >= omega1);
  
  ier = 1;
  phi = geogr_to_geoce(phi_g);
  ier = sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
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
      a = I0 * eccentricity * DAY_LENGTH / (2.0 * Pi);
      b1 = sin(phi) * sin(delta) * (omega2 - omega1);
      b2 = cos(phi) * cos(delta) * (sin(omega2) - sin(omega1));
      c = a * (b1 + b2);
      if (c < 0.0)
        *G0_12 = 0.0;
      else
        *G0_12 = c;
    }

  return (ier);
}

int G0_day(double phi_g, double eccentricity, double delta, 
  double *G0d)
{
  int ier;
  double omega_sr, omega_ss, a, b;
  double phi;

  ier = 1;
  phi = geogr_to_geoce(phi_g);
  ier = sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = I0 * eccentricity * DAY_LENGTH / Pi;
  /*
   * b = cos(phi) * cos(delta) * (sin(omega_ss) - omega_ss * cos(omega_ss));
   */
  b = sin(phi) * sin(delta) * omega_ss + cos(phi) * cos(delta) * sin(omega_ss);
  *G0d = a * b;

  return (ier);
}

int G0_hours_profile(double phi_g, double eccentricity, double delta, 
  double *G0h)
{
  int ier, i;
  double omega_sr, omega_ss, a, b1, b2;
  double phi;
  double t1, t2, omega1, omega2;

  ier = 1;
  phi = geogr_to_geoce(phi_g);
  ier = sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = I0 * eccentricity * DAY_LENGTH / (2.0 * Pi);
  b1 = sin(phi) * sin(delta);
  b2 = cos(phi) * cos(delta);

  for (i = 0; i < 24; i++)
    {
      t1 = (double) (i + 1) - 1.0;
      omega1 = solar_hour_angle(t1);
      t2 = (double) (i + 1);
      omega2 = solar_hour_angle(t2);

      if ((omega2 < omega_sr) || (omega1 > omega_ss))
        G0h[i] = 0.0;
      else
        {
          if (omega1 < omega_sr)
            omega1 = omega_sr;
          if (omega2 > omega_ss)
            omega2 = omega_ss;
          G0h[i] = a * (b1 * (omega2 - omega1) + b2 * (sin(omega2)
              - sin(omega1)));
        }
    }

  return (ier);
}

int G0_hour(double phi_g, double eccentricity, double delta, double t, 
  double *G0h)
{
  int ier;
  double omega_sr, omega_ss, a, b1, b2;
  double t1, t2, omega1, omega2;
  double phi;

  ier = 1;
  phi = geogr_to_geoce(phi_g);
  ier = sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
  if (ier != 0)
    return (ier);
  a = I0 * eccentricity * DAY_LENGTH / (2.0 * Pi);
  b1 = sin(phi) * sin(delta);
  b2 = cos(phi) * cos(delta);

  t1 = t - 1.0;
  omega1 = solar_hour_angle(t1);
  t2 = t;
  omega2 = solar_hour_angle(t2);

  if (omega1 < omega_sr)
    omega1 = omega_sr;
  if (omega2 < omega_sr)
    omega2 = omega_sr;

  if (omega2 <= omega1)
    *G0h = 0.0;
  else
    {
      *G0h = a * (b1 * (omega2 - omega1) + b2 * (sin(omega2) - sin(omega1)));
      if (*G0h < 0.0)
        *G0h = 0.0;
    }

  return (ier);
}

/***********************************************/
/* MONTHLY AVERAGES OF SOLAR INPUTS PARAMETERS */
/***********************************************/
int monthly_averages(int month_number, int year_number, 
        double phi_g, double lambda, double gamma_riset, 
        double *day_angle_m, double *delta_m, double *omega_ss_m,
        double *S0_m, double *eccentricity_m,
        double *G0d_m, double *G0h_m)
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

  ier = 0;
  number_days_month = nbdays_month(year_number, month_number);

  for (day_of_month = 1; day_of_month <= number_days_month; day_of_month++)
    {
      ier = make_julian_day(day_of_month, month_number, year_number,
          &julian_day);
      if (ier == 0)
        day_angle = Day_Angle(julian_day);
      if (ier == 0)
          delta = declination_sun(year_number, julian_day, lambda);
      if (ier == 0)
        ier = sunrise_hour_angle(phi_g, delta, gamma_riset, &omega_sr,
            &omega_ss);
      if (ier == 0)
        ier = timerise_daylength(omega_sr, omega_ss, &t_sr, &t_ss, &S0);
      if (ier == 0)
        eccentricity = corr_distance(day_angle);
      if (ier == 0)
        ier = G0_day(phi_g, eccentricity, delta, &G0d);
      if (ier == 0)
        ier = G0_hours_profile(phi_g, eccentricity, delta, G0h);
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
/* YEARLY AVERAGES OF MONTHLY SOLAR PARAMETERS      */
/* (LONG TERM MEANS OF MONTHLY MEANS OF DAILY SUMS) */
/****************************************************/
int yearly_averages(int month_number, int year_start, int year_end, 
        double phi_g, double lambda, double gamma_riset, 
        double *day_angle_y, double *delta_y, double *omega_ss_y, 
        double *S0_y, double *eccentricity_y, 
        double *G0d_y, double *G0h_y)
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
      ier = monthly_averages(month_number, year_number, phi_g, lambda,
          gamma_riset, &day_angle_m, &delta_m, &omega_ss_m, &S0_m,
          &eccentricity_m, &G0d_m, G0h_m);
      if (ier != 0)
        return (ier);

      *day_angle_y = *day_angle_y + day_angle_m;
      *delta_y = *delta_y + delta_m;
      *omega_ss_y = *omega_ss_y + omega_ss_m;
      *S0_y = *S0_y + S0_m;
      *eccentricity_y = *eccentricity_y + eccentricity_m;
      *G0d_y = *G0d_y + G0d_m;
      for (i = 1; i <= 24; i++)
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
/* SOLAR INPUTS PARAMETERS FOR A CERTAIN DAY */
/*********************************************/
int solar_parameters_day(int day_of_month, int month_number, int year_number,
        double phi_g, double lambda, double gamma_riset, 
        double *day_angle, double *delta, double *omega_ss, 
        double *S0, double *eccentricity, 
        double *G0d, double *G0h)
{
  int ier, julian_day;
  double omega_sr, t_sr, t_ss;

  ier = 1;
  ier = make_julian_day(day_of_month, month_number, year_number, &julian_day);
  if (ier == 0)
    *day_angle = Day_Angle(julian_day);
  if (ier == 0)
    *delta = declination_sun(year_number, julian_day, lambda);
  if (ier == 0)
    ier = sunrise_hour_angle(phi_g, *delta, gamma_riset, &omega_sr, omega_ss);
  if (ier == 0)
    ier = timerise_daylength(omega_sr, *omega_ss, &t_sr, &t_ss, S0);
  if (ier == 0)
    *eccentricity = corr_distance(*day_angle);
  if (ier == 0)
    ier = G0_day(phi_g, *eccentricity, *delta, G0d);
  if (ier == 0 && *G0d > 0.0)
    ier = G0_hours_profile(phi_g, *eccentricity, *delta, G0h);

  return (ier);
}

/******************************************************************/
/* SOLAR INPUTS PARAMETERS FOR MONTHLY AVERAGE IRRADIATION MODELS */
/******************************************************************/
int solar_parameters_avg(int month_number, 
			   double phi_g, double gamma_riset,
			   double *day_angle_avg,
			   double *delta_avg,
			   double *omega_ss_avg, double *S0_avg,
			   double *eccentricity_avg,
			   double *G0d_avg, double *G0h_avg)
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
  type_use = 0; /* for estimating monthly mean global solar radiation */

  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 0)
        julian_day = tab_julian_day[month_number - 1];
      ier = 0;
    }
  if (ier == 0)
    *day_angle_avg = Day_Angle(julian_day);
  if (ier == 0)
    ier = declination_sun_month(month_number, type_use, delta_avg);
  if (ier == 0)
    ier = sunrise_hour_angle(phi_g, *delta_avg, gamma_riset, &omega_sr,
        omega_ss_avg);
  if (ier == 0)
    ier = timerise_daylength(omega_sr, *omega_ss_avg, &t_sr, &t_ss, S0_avg);
  if (ier == 0)
    *eccentricity_avg = corr_distance(*day_angle_avg);
  if (ier == 0)
    ier = G0_day(phi_g, *eccentricity_avg, *delta_avg, G0d_avg);
  if (ier == 0)
    ier = G0_hours_profile(phi_g, *eccentricity_avg, *delta_avg, G0h_avg);

  return (ier);
}

/******************************************************************/
/* SOLAR INPUTS PARAMETERS FOR MONTHLY MAXIMUM IRRADIATION MODELS */
/******************************************************************/
int solar_parameters_max(int month_number, 
			   double phi_g, double gamma_riset,
			   double *day_angle_max,
			   double *delta_max,
			   double *omega_ss_max, double *S0_max,
			   double *eccentricity_max,
			   double *G0d_max, double *G0h_max)
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
  type_use = 1; /* for estimating monthly mean global solar radiation */

  if ((type_use >= 0) && (type_use < 2))
    {
      if (type_use == 1)
        julian_day = tab_julian_day_max[month_number - 1];
      ier = 0;
    }
  if (ier == 0)
    *day_angle_max = Day_Angle(julian_day);
  if (ier == 0)
    ier = declination_sun_month(month_number, type_use, delta_max);
  if (ier == 0)
    ier = sunrise_hour_angle(phi_g, *delta_max, gamma_riset, &omega_sr,
        omega_ss_max);
  if (ier == 0)
    ier = timerise_daylength(omega_sr, *omega_ss_max, &t_sr, &t_ss, S0_max);
  if (ier == 0)
    *eccentricity_max = corr_distance(*day_angle_max);
  if (ier == 0)
    ier = G0_day(phi_g, *eccentricity_max, *delta_max, G0d_max);
  if (ier == 0)
    ier = G0_hours_profile(phi_g, *eccentricity_max, *delta_max, G0h_max);

  return (ier);
}

int intervals_omega_tilted_plane(double phi_g, double delta, double omega_ss, double beta, double alpha, 
    double *v_om, int *p_nb)
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

  v_om[0] = -999.0;
  v_om[1] = -999.0;
  v_om[2] = -999.0;
  v_om[3] = -999.0;
  *p_nb = 0;

  omega_sr = -omega_ss;

  /* Le code ci-dessous pour optimiser les plages de calcul effectif entre omega_sr et omega_ss
   * semble marcher pour l'hémisphère Nord mais pas pour l'hémisphère Sud : en attendant, on shunte
   * l'optimisation en forçant le calcul complet entre le lever et le coucher du soleil
   *  */
  v_om[0] = omega_sr;
  v_om[1] = omega_ss;
  *p_nb = 2;
  return (ierr);

  /* Fin du shunte */

  if (fabs(omega_ss) < precision)
    {
      return ierr;
    }

  phi = geogr_to_geoce(phi_g);

  A = cos(delta) * (cos(phi) * cos(beta) + sin(phi) * sin(beta) * cos(alpha));

  B = cos(delta) * sin(beta) * sin(alpha);
  C = sin(delta) * (sin(phi) * cos(beta) - cos(phi) * sin(beta) * cos(alpha));

  nzA = (int) (fabs(A) > epsilon_z);
  nzB = (int) (fabs(B) > epsilon_z);
  nzC = (int) (fabs(C) > epsilon_z);

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
    sAB = sqrt(mAB);

    d = fabs(C) / sAB;
    if (fabs(d - 1) < precision)
      {

        if (C > 0.0)
          {
            wa = -asin(B);
          }
        else
          {
            wa = asin(B);
          }

        if ((wa < omega_sr + precision) | (wa > omega_ss - precision))
          {
            wt = (omega_sr + omega_ss) / 2.0;
            cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
            if (cosTheta_wt > precision)
              {
                /*
                 * ANGLES = [omega_sr omega_ss];
                 */
                v_om[0] = omega_sr;
                v_om[1] = omega_ss;
              }
            else if (fabs(cosTheta_wt) < precision)
              {
                ierr = 1;
              }
          }
        else
          {
            /*
             * Solution soit [omega_sr wa] ou [wa omega_ss] définir Test au
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
            cosTheta_wm = cos(wm) * A + sin(wm) * B + C;
            cosTheta_wt = cos(wt) * A + sin(wt) * B + C;

            if (fabs(cosTheta_wm) > fabs(cosTheta_wt))
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
        cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
        if (cosTheta_wt > precision)
          {
            /*
             * ANGLES = [omega_sr omega_ss];
             */
            v_om[0] = omega_sr;
            v_om[1] = omega_ss;
          }
        else if (fabs(cosTheta_wt) < precision)
          {
            ierr = 1;
          }
      }
    else
      {
        if (cas != 6)
          {

            r1 = fabs(A) * sqrt(mAB - C * C);
            sin_w1 = -(B * C + r1) / mAB;
            sin_w2 = -(B * C - r1) / mAB;

            r2 = fabs(B) * sqrt(mAB - C * C);
            cos_w1 = -(A * C + r2) / mAB;
            cos_w2 = -(A * C - r2) / mAB;

            if (fabs(sin_w1 * sin_w1 + cos_w1 * cos_w1 - 1.0) > epsilon_z)
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

            wa = atan2(sin_wa, cos_wa);
            wb = atan2(sin_wb, cos_wb);
          }
        else
          {
            wa = atan(-A / B);
            wb = Pi + wa;
            if (wb > Pi)
              {
                wb = wb - 2.0 * Pi;
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
            cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
            if (cosTheta_wt > precision)
              {
                /*
                 * ANGLES = [omega_sr omega_ss];
                 */
                v_om[0] = omega_sr;
                v_om[1] = omega_ss;
              }
            else if (fabs(cosTheta_wt) < precision)
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
            cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
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
            cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
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
            cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
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

  return (ierr);

}
