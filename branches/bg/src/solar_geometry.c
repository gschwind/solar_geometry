/*-------------------------------------------------------------------------*/
/*              GEOMETRY OF SOLAR BEAM, A NEAR POINTSOURCE                 */
/*-------------------------------------------------------------------------*/
/*                        ECOLE DES MINES DE PARIS                         */
/*        CENTRE D'ENERGETIQUE - GROUPE TELEDETECTION & MODELISATION       */
/*                       Rue Claude Daunesse, BP 207                       */
/*                   06904 Sophia Antipolis cedex, FRANCE                  */
/*          Tel (+33) 04 93 95 74 49     Fax (+33) 04 93 95 75 35          */
/*                   E-mail : lucien.wald@mines-paristech.fr               */
/*-------------------------------------------------------------------------*/
/*   L. Wald - O. Bauer - February 1997                                    */
/*   modified 8 July 2004 L. Wald for geocentric - geographic lat          */
/*-------------------------------------------------------------------------*/

#define __C_solar_geometry

#include "solar_geometry.h"

#include <math.h>
#include <stdio.h>

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define IS_LEAP_YEAR(y) (((((y) % 4) == 0) && (((y) % 100) != 0)) || (((y) % 400) == 0))

/*
 * 0. BASIC PARAMETERS
 * --------------------
 */

int ymd_to_day_of_year(int year, int month, int day_of_month) {
	static int const tab[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
			304, 334 };
	int day_of_year = day_of_month + tab[month - 1];
	if (IS_LEAP_YEAR(year) && (month > 2))
		day_of_year = day_of_year + 1;
	return day_of_year;
}

void day_of_year_to_ymd(int year, int day_of_year, int * day_of_month,
		int * month) {
	static int const tab0[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
			304, 334 };
	static int const tab1[12] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274,
			305, 335 };
	int const * tab;
	int m;
	/* leap year */
	if (IS_LEAP_YEAR(year)) {
		tab = tab1;
	} else {
		tab = tab0;
	}

	for (m = 11; m >= 0; ++m) {
		if (day_of_year > tab[m]) {
			*month = m + 1;
			*day_of_month = day_of_year - tab[m];
			break;
		}
	}
}

int nbdays_month(int year_number, int month_number) {
	static int tab_nbdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30,
			31 };
	if (IS_LEAP_YEAR(year_number) && (month_number == 2))
		return tab_nbdays[month_number - 1] + 1;
	else
		return tab_nbdays[month_number - 1];
}

char const * const get_month_name(int month) {
	char const * const tab_name[13] = { "", "jan", "feb", "mar", "apr", "may",
			"jun", "jul", "aug", "sep", "oct", "nov", "dec" };
	return tab_name[month];
}

double get_day_angle(int day_of_year) {
	return ((double) day_of_year) * 2.0 * M_PI / 365.2422;
}

double declination_sun(int year, int day_of_year, double lambda) {
	double const b1 = 0.0064979;
	double const b2 = 0.4059059;
	double const b3 = 0.0020054;
	double const b4 = -0.0029880;
	double const b5 = -0.0132296;
	double const b6 = 0.0063809;
	double const b7 = 0.0003508;

	double w0, n0, t1, wt;

	/*
	 * n0 : spring-equinox time expressed in days from the beginning of the year i.e.
	 * the time in decimal days elapsing from 00:00 hours Jan 1st to the spring equinox
	 * at Greenwich in a given year
	 *
	 * t1 : time in days, from the spring equinox
	 *
	 * 0.5 represents the decimal day number at noon on Jan 1st at Greenwich
	 */
	n0 = 78.8946 + 0.2422 * (year - 1957) - (int) (0.25 * (year - 1957));
	t1 = -0.5 - lambda / (2 * M_PI) - n0;
	w0 = 2 * M_PI / 365.2422;
	wt = w0 * (day_of_year + t1);
	return b1 + b2 * sin(wt) + b3 * sin(2 * wt) + b4 * sin(3 * wt) + b5 * cos(
			wt) + b6 * cos(2 * wt) + b7 * cos(3 * wt);
}

double declination_sun_month(int month, int type_use) {
	double const deg_rad = (M_PI / 180.0); /* converts decimal degrees into radians */
	int tab_day_of_year[12] = { 17, 46, 75, 105, 135, 162, 198, 228, 259, 289,
			319, 345 };
	int tab_day_of_year_max[12] = { 29, 57, 89, 119, 150, 173, 186, 217, 248,
			278, 309, 339 };
	int day_of_year;
	double day_angle;

	if ((type_use >= 0) && (type_use < 2)) {
		if (type_use == 0)
			day_of_year = tab_day_of_year[month - 1];
		if (type_use == 1)
			day_of_year = tab_day_of_year_max[month - 1];
	}

	day_angle = get_day_angle(day_of_year);
	double const c1 = 0.3978;
	double const c2 = 80.2 * deg_rad; /* 1.4000 in SSA manual */
	double const c3 = 1.92 * deg_rad; /* 0.0355 in SSA manual */
	double const c4 = 2.80 * deg_rad; /* 0.0489 in SSA manual */

	return asin(c1 * sin(day_angle - c2 + c3 * sin(day_angle - c4)));
}

double solar_hour_angle(double t) {
	return (t - 12.0) * M_PI / 12.0;
}

double omega_to_LAT(double omega) {
	return 12.0 * (1.0 + omega / M_PI);
}

double geogr_to_geoce(double phi_g) {
	double const CC = 0.99330552; /* Correction factor for converting geographic */
	/*  into geocentric latitude. CC=(Rpole/Requator)**2
	 * Rpole=6356.752, Requator=6378.137
	 */
	if ((phi_g >= -(M_PI / 2.0 - 0.0002)) || (phi_g <= (M_PI / 2.0 - 0.0002)))
		return atan(tan(phi_g) * CC);
	else
		return phi_g;
}

double solar_hour_angle_h(double phi_g, double delta, double t) {
	double omega_sr, omega_ss, omega1, omega2;

	sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
	omega1 = (t - 1.0 - 12.0) * M_PI / 12.0;
	if (omega1 < omega_sr)
		omega1 = omega_sr;

	omega2 = (t - 12.0) * M_PI / 12.0;
	if (omega2 > omega_ss)
		omega2 = omega_ss;
	return (omega1 + omega2) / 2.0;
}

/*
 * 1. SUNRISE, SUNSET AND DAYLENGTH
 * ---------------------------------
 */

void sunrise_hour_angle(double phi_g, double delta, double gamma_riset,
		double *omega_sr, double *omega_ss) {
	static double const deg_rad = (M_PI / 180.0); /* converts decimal degrees into radians */
	double horizon, max_delta, cos_omegas, omegas;
	double phi;

	horizon = (-50.0 / 60.0) * deg_rad; /* horizon, -50' in radians */
	if (gamma_riset >= horizon)
		horizon = gamma_riset;

	phi = geogr_to_geoce(phi_g);
	max_delta = 23.45 * deg_rad;

	cos_omegas = (sin(horizon) - (sin(phi) * sin(delta))) / (cos(phi) * cos(
			delta));

	if (fabs(cos_omegas) < 1.0)
		omegas = acos(cos_omegas);
	if (cos_omegas >= 1.0) /* the sun is always below the horizon : polar night */
		omegas = 0.0;
	if (cos_omegas <= -1.0) /* the sun is always above the horizon : polar day */
		omegas = M_PI;

	*omega_sr = -omegas;
	*omega_ss = omegas;
}

void timerise_daylength(double omega_sr, double omega_ss, double *t_sr,
		double *t_ss, double *S0) {
	/*
	 * alternative way
	 *
	 * ier = omega_to_LAT(omega_sr,&t_sr); if(ier == 0) ier ==
	 * omega_to_LAT(omega_ss,&t_ss); if(ier != 0) return(ier);
	 */
	*t_sr = 12.0 + omega_sr * 12.0 / M_PI;
	*t_ss = 12.0 + omega_ss * 12.0 / M_PI;
	*S0 = *t_ss - *t_sr;
}

/*
 * 2. CHANGING THE TIME SYSTEM
 * ----------------------------
 */

double LMT_to_LAT(double day_angle, double lambda, double lambda_ref,
		int summer_corr) {
	double const deg_rad = (M_PI / 180.0); /* converts decimal degrees into radians */
	double const a1 = -0.128;
	double const a2 = -0.165;
	double const a3 = 2.80 * deg_rad;
	double const a4 = 19.70 * deg_rad;
	double ET;
	ET = a1 * sin(day_angle - a3) + a2 * sin(2.0 * day_angle + a4);
	return ET + ((lambda - lambda_ref) * 12.0 / M_PI) - (double) summer_corr;
}

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
double UT_to_LAT(double UT, double day_angle, double lambda) {
	const double deg_rad = (M_PI / 180.0); /* converts decimal degrees into radians */
	double const a1 = -0.128;
	double const a2 = -0.165;
	double const a3 = 2.80 * deg_rad;
	double const a4 = 19.70 * deg_rad;
	double ET, LAT;

	ET = a1 * sin(day_angle - a3) + a2 * sin(2.0 * day_angle + a4);
	LAT = UT + ET + (lambda * 12.0 / M_PI);
	if (LAT < 0)
		LAT += 24.0;
	if (LAT > 24.0)
		LAT -= 24.0;
	return LAT;
}

/*
 * 3. POSITION OF THE SUN IN THE SKY
 * ----------------------------------
 */

void elevation_zenith_sun(double phi_g, double delta, double omega,
		double *gamma, double *theta) {
	double omega_sr, omega_ss;
	double phi;
	phi = geogr_to_geoce(phi_g);
	sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
	if ((omega < omega_sr) || (omega > omega_ss))
		*gamma = 0.0;
	else
		*gamma = asin(sin(phi) * sin(delta) + cos(phi) * cos(delta)
				* cos(omega));
	if (*gamma < 0.0)
		*gamma = 0.0;

	*theta = (M_PI / 2.0) - *gamma;
}

double azimuth_sun(double phi_g, double delta, double omega, double gamma) {
	double cos_as, sin_as, x;
	double phi;
	phi = geogr_to_geoce(phi_g);
	cos_as = (sin(phi) * sin(gamma) - sin(delta)) / (cos(phi) * cos(gamma));
	if (phi < 0.0)
		cos_as = -cos_as; /* Southern hemisphere */
	sin_as = cos(delta) * sin(omega) / cos(gamma);

	if (cos_as > 1.0)
		cos_as = 1.0;

	if (cos_as < -1.0)
		cos_as = -1.0;

	x = acos(cos_as);
	if (sin_as >= 0.0)
		return x;
	else
		return -x;
}

/*
 * 4. EXTRATERRESTRIAL IRRADIATION
 * --------------------------------
 */

double corr_distance(double day_angle) {
	double const deg_rad = (M_PI / 180.0); /* converts decimal degrees into radians */
	double const a = 2.80 * deg_rad;
	return 1.0 + 0.03344 * cos(day_angle - a);
}

double G0_normal(double I0j, double theta) {
	return I0j * cos(theta);
}

double G0_general(double phi_g, double eccentricity, double delta,
		double omega1, double omega2) {
	double omega_sr, omega_ss, a, b1, b2, c;
	double phi;

	phi = geogr_to_geoce(phi_g);
	sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);

	if (omega1 < omega_sr)
		omega1 = omega_sr;
	if (omega2 < omega_sr)
		omega2 = omega_sr;
	if (omega1 > omega_ss)
		omega1 = omega_ss;
	if (omega2 > omega_ss)
		omega2 = omega_ss;

	if (omega2 <= omega1) {
		return 0.0;
	} else {
		a = I0 * eccentricity * DAY_LENGTH / (2.0 * M_PI);
		b1 = sin(phi) * sin(delta) * (omega2 - omega1);
		b2 = cos(phi) * cos(delta) * (sin(omega2) - sin(omega1));
		c = a * (b1 + b2);
		if (c < 0.0)
			return 0.0;
		else
			return c;
	}

}

double G0_day(double phi_g, double eccentricity, double delta) {
	double omega_sr, omega_ss, a, b;
	double phi;
	phi = geogr_to_geoce(phi_g);
	sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
	a = I0 * eccentricity * DAY_LENGTH / M_PI;
	/*
	 * b = cos(phi) * cos(delta) * (sin(omega_ss) - omega_ss * cos(omega_ss));
	 */
	b = sin(phi) * sin(delta) * omega_ss + cos(phi) * cos(delta)
			* sin(omega_ss);
	return a * b;
}

void G0_hours_profile(double phi_g, double eccentricity, double delta,
		double * G0h) {
	int i;
	double omega_sr, omega_ss, a, b1, b2;
	double phi;
	double t1, t2, omega1, omega2;

	phi = geogr_to_geoce(phi_g);
	sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
	a = I0 * eccentricity * DAY_LENGTH / (2.0 * M_PI);
	b1 = sin(phi) * sin(delta);
	b2 = cos(phi) * cos(delta);

	for (i = 0; i < 24; i++) {
		t1 = (double) (i + 1) - 1.0;
		omega1 = solar_hour_angle(t1);
		t2 = (double) (i + 1);
		omega2 = solar_hour_angle(t2);
		if ((omega2 < omega_sr) || (omega1 > omega_ss))
			G0h[i] = 0.0;
		else {
			if (omega1 < omega_sr)
				omega1 = omega_sr;
			if (omega2 > omega_ss)
				omega2 = omega_ss;
			G0h[i] = a * (b1 * (omega2 - omega1) + b2 * (sin(omega2) - sin(
					omega1)));
		}
	}
}

double G0_hour(double phi_g, double eccentricity, double delta, double t) {
	double omega_sr, omega_ss, a, b1, b2;
	double t1, t2, omega1, omega2;
	double phi, G0h;

	phi = geogr_to_geoce(phi_g);
	sunrise_hour_angle(phi_g, delta, 0.0, &omega_sr, &omega_ss);
	a = I0 * eccentricity * DAY_LENGTH / (2.0 * M_PI);
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
		G0h = 0.0;
	else {
		G0h = a * (b1 * (omega2 - omega1) + b2 * (sin(omega2) - sin(omega1)));
		if (G0h < 0.0)
			G0h = 0.0;
	}
	return G0h;
}

/*
 * 4.1. MONTHLY AVERAGES OF SOLAR INPUTS PARAMETERS
 * -------------------------------------------------
 */

void monthly_averages(int month, int year, double phi_g, double lambda,
		double gamma_riset, double *day_angle_m, double *delta_m,
		double *omega_ss_m, double *S0_m, double *eccentricity_m,
		double *G0d_m, double *G0h_m) {
	int i, day_of_month, number_days_month, julian_day;
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

	number_days_month = nbdays_month(year, month);

	for (day_of_month = 1; day_of_month <= number_days_month; day_of_month++) {
		julian_day = ymd_to_day_of_year(day_of_month, month, year);
		day_angle = get_day_angle(julian_day);
		delta = declination_sun(year, julian_day, lambda);
		sunrise_hour_angle(phi_g, delta, gamma_riset, &omega_sr, &omega_ss);
		timerise_daylength(omega_sr, omega_ss, &t_sr, &t_ss, &S0);
		eccentricity = corr_distance(day_angle);
		G0d = G0_day(phi_g, eccentricity, delta);
		G0_hours_profile(phi_g, eccentricity, delta, G0h);

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
}

/*
 * 4.2. YEARLY AVERAGES OF MONTHLY SOLAR PARAMETERS
 * (LONG TERM MEANS OF MONTHLY MEANS OF DAILY SUMS)
 * -------------------------------------------------
 */

void yearly_averages(int month, int year_start, int year_end, double phi_g,
		double lambda, double gamma_riset, double *day_angle_y,
		double *delta_y, double *omega_ss_y, double *S0_y,
		double *eccentricity_y, double *G0d_y, double *G0h_y) {
	int i, year;
	double day_angle_m, delta_m, omega_ss_m, S0_m, eccentricity_m, G0d_m,
			G0h_m[24];
	double number_of_years;

	number_of_years = (double) (year_end - year_start + 1.);

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

	for (year = year_start; year <= year_end; year++) {
		monthly_averages(month, year, phi_g, lambda, gamma_riset, &day_angle_m,
				&delta_m, &omega_ss_m, &S0_m, &eccentricity_m, &G0d_m, G0h_m);

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
}

/*
 * 4.3. SOLAR INPUTS PARAMETERS FOR A CERTAIN DAY
 * -----------------------------------------------
 */

void solar_parameters_day(int day_of_month, int month_number, int year_number,
		double phi_g, double lambda, double gamma_riset, double *day_angle,
		double *delta, double *omega_ss, double *S0, double *eccentricity,
		double *G0d, double *G0h) {
	int day_of_year;
	double omega_sr, t_sr, t_ss;

	day_of_year = ymd_to_day_of_year(day_of_month, month_number, year_number);
	*day_angle = get_day_angle(day_of_year);
	*delta = declination_sun(year_number, day_of_year, lambda);
	sunrise_hour_angle(phi_g, *delta, gamma_riset, &omega_sr, omega_ss);
	timerise_daylength(omega_sr, *omega_ss, &t_sr, &t_ss, S0);
	*eccentricity = corr_distance(*day_angle);
	*G0d = G0_day(phi_g, *eccentricity, *delta);
	if (*G0d > 0.0)
		G0_hours_profile(phi_g, *eccentricity, *delta, G0h);

}

/*
 * 4.4. SOLAR INPUTS PARAMETERS FOR MONTHLY AVERAGE IRRADIATION MODELS
 * --------------------------------------------------------------------
 */

void solar_parameters_avg(int month, double phi_g, double gamma_riset,
		double *day_angle_avg, double *delta_avg, double *omega_ss_avg,
		double *S0_avg, double *eccentricity_avg, double *G0d_avg,
		double *G0h_avg) {
	int day_of_year;
	double omega_sr, t_sr, t_ss;

	/*
	 * recommended values of day number for estimating monthly mean global solar
	 * radiation
	 */
	int tab_day_of_year[12] = { 17, 46, 75, 105, 135, 162, 198, 228, 259, 289,
			319, 345 };
	day_of_year = tab_day_of_year[month - 1];
	*day_angle_avg = get_day_angle(day_of_year);
	*delta_avg = declination_sun_month(month, 0);
	sunrise_hour_angle(phi_g, *delta_avg, gamma_riset, &omega_sr, omega_ss_avg);
	timerise_daylength(omega_sr, *omega_ss_avg, &t_sr, &t_ss, S0_avg);
	*eccentricity_avg = corr_distance(*day_angle_avg);
	*G0d_avg = G0_day(phi_g, *eccentricity_avg, *delta_avg);
	G0_hours_profile(phi_g, *eccentricity_avg, *delta_avg, G0h_avg);
}

/*
 * 4.5. SOLAR INPUTS PARAMETERS FOR MONTHLY MAXIMUM IRRADIATION MODELS
 * --------------------------------------------------------------------
 */

void solar_parameters_max(int month, double phi_g, double gamma_riset,
		double *day_angle_max, double *delta_max, double *omega_ss_max,
		double *S0_max, double *eccentricity_max, double *G0d_max,
		double *G0h_max) {
	int ier, day_of_year;
	double omega_sr, t_sr, t_ss;

	/*
	 * recommended values of day number for estimating monthly mean maximum global solar
	 * radiation
	 */
	int tab_day_of_year_max[12] = { 29, 57, 89, 119, 150, 173, 186, 217, 248,
			278, 309, 339 };

	day_of_year = tab_day_of_year_max[month - 1];
	*day_angle_max = get_day_angle(day_of_year);
	*delta_max = declination_sun_month(month, 1);
	sunrise_hour_angle(phi_g, *delta_max, gamma_riset, &omega_sr, omega_ss_max);
	timerise_daylength(omega_sr, *omega_ss_max, &t_sr, &t_ss, S0_max);
	*eccentricity_max = corr_distance(*day_angle_max);
	*G0d_max = G0_day(phi_g, *eccentricity_max, *delta_max);
	G0_hours_profile(phi_g, *eccentricity_max, *delta_max, G0h_max);
}


void intervals_omega_tilted_plane(double phi_g, double delta, double omega_ss,
		double beta, double alpha, double *v_om, int *p_nb) {
	double A, B, C;
	double const precision = 1e-4;
	double const epsilon_z = 1e-6;
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

	phi = geogr_to_geoce(phi_g);

	v_om[0] = -999.0;
	v_om[1] = -999.0;
	v_om[2] = -999.0;
	v_om[3] = -999.0;
	*p_nb = 0;

	if (fabs(omega_ss) < precision) {
		return;
	}

	omega_sr = -omega_ss;

	A = cos(delta) * (cos(phi) * cos(beta) + sin(phi) * sin(beta) * cos(alpha));

	B = cos(delta) * sin(beta) * sin(alpha);
	C = sin(delta) * (sin(phi) * cos(beta) - cos(phi) * sin(beta) * cos(alpha));

	nzA = (int) (fabs(A) > epsilon_z);
	nzB = (int) (fabs(B) > epsilon_z);
	nzC = (int) (fabs(C) > epsilon_z);

	cas = nzC + nzB * 2 + nzA * 4;

	switch (cas) {
	case 0:
		/*
		 * ANGLES = [omega_sr omega_ss];
		 */
		v_om[0] = omega_sr;
		v_om[1] = omega_ss;
		break;
	case 1:
		if (C > 0) {
			/*
			 * ANGLES = [omega_sr omega_ss];
			 */
			v_om[0] = omega_sr;
			v_om[1] = omega_ss;
		}
		break;
	case 2:
		if (B > 0) {
			/*
			 * ANGLES = [0 omega_ss];
			 */
			v_om[0] = 0.0;
			v_om[1] = omega_ss;
		} else {
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
		if (fabs(d - 1) < precision) {

			if (C > 0.0) {
				wa = -asin(B);
			} else {
				wa = asin(B);
			}

			if ((wa < omega_sr + precision) | (wa > omega_ss - precision)) {
				wt = (omega_sr + omega_ss) / 2.0;
				cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
				if (cosTheta_wt > precision) {
					/*
					 * ANGLES = [omega_sr omega_ss];
					 */
					v_om[0] = omega_sr;
					v_om[1] = omega_ss;
				} else if (fabs(cosTheta_wt) < precision) {

				}
			} else {
				/*
				 * Solution soit [omega_sr wa] ou [wa omega_ss] ���d���inir Test au
				 * centre et au bord omega_sr ou omega_ss du plus grand intervalle
				 */
				if (omega_ss - wa > wa - omega_sr) {
					wt = (omega_ss + wa) / 2;
					wm = omega_ss;
				} else {
					wt = (omega_sr + wa) / 2;
					wm = omega_sr;
				}
				cosTheta_wm = cos(wm) * A + sin(wm) * B + C;
				cosTheta_wt = cos(wt) * A + sin(wt) * B + C;

				if (fabs(cosTheta_wm) > fabs(cosTheta_wt)) {
					cosTheta_wt = cosTheta_wm;
				}
				if (cosTheta_wt > precision) {
					/*
					 * ANGLES = [wa omega_ss];
					 */
					v_om[0] = wa;
					v_om[1] = omega_ss;
				} else if (cosTheta_wt < -precision) {
					/*
					 * ANGLES = [omega_sr wa];
					 */
					v_om[0] = omega_sr;
					v_om[1] = wa;
				} else {

				}
			}
		} else if (d >= 1 + precision) {

			wt = (omega_sr + omega_ss) / 2.0;
			cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
			if (cosTheta_wt > precision) {
				/*
				 * ANGLES = [omega_sr omega_ss];
				 */
				v_om[0] = omega_sr;
				v_om[1] = omega_ss;
			} else if (fabs(cosTheta_wt) < precision) {

			}
		} else {
			if (cas != 6) {

				r1 = fabs(A) * sqrt(mAB - C * C);
				sin_w1 = -(B * C + r1) / mAB;
				sin_w2 = -(B * C - r1) / mAB;

				r2 = fabs(B) * sqrt(mAB - C * C);
				cos_w1 = -(A * C + r2) / mAB;
				cos_w2 = -(A * C - r2) / mAB;

				if (fabs(sin_w1 * sin_w1 + cos_w1 * cos_w1 - 1.0) > epsilon_z) {
					sin_wa = sin_w1;
					cos_wa = cos_w2;
					sin_wb = sin_w2;
					cos_wb = cos_w1;
				} else {
					sin_wa = sin_w1;
					cos_wa = cos_w1;
					sin_wb = sin_w2;
					cos_wb = cos_w2;
				}

				wa = atan2(sin_wa, cos_wa);
				wb = atan2(sin_wb, cos_wb);
			} else {
				wa = atan(-A / B);
				wb = M_PI + wa;
				if (wb > M_PI) {
					wb = wb - 2.0 * M_PI;
				}
			}

			wmin = MIN (wa, wb);
			wmax = MAX (wa, wb);
			wa = wmin;
			wb = wmax;

			pa = 0;
			if (wa < omega_sr + precision) {
				pa = -1;
			} else if (wa > omega_ss - precision) {
				pa = 1;
			}

			pb = 0;
			if (wb < omega_sr + precision) {
				pb = -1;
			} else if (wb > omega_ss - precision) {
				pb = 1;
			}

			if ((pa != 0) && (pb != 0)) {
				/*
				 * [wa omega_sr omega_ss wb]
				 */
				wt = (omega_sr + omega_ss) / 2.0;
				cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
				if (cosTheta_wt > precision) {
					/*
					 * ANGLES = [omega_sr omega_ss];
					 */
					v_om[0] = omega_sr;
					v_om[1] = omega_ss;
				} else if (fabs(cosTheta_wt) < precision) {

				}
			} else if ((pa == 0) && (pb == 0)) {
				/*
				 * [omega_sr wa wb omega_ss]
				 */
				wt = (wa + wb) / 2.0;
				cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
				if (cosTheta_wt > precision) {
					/*
					 * ANGLES = [wa wb];
					 */
					v_om[0] = wa;
					v_om[1] = wb;
				} else if (cosTheta_wt < -precision) {
					/*
					 * ANGLES = [omega_sr wa wb omega_ss];
					 */
					v_om[0] = omega_sr;
					v_om[1] = wa;
					v_om[2] = wb;
					v_om[3] = omega_ss;
				} else {

				}
			} else if (pa == 0) {
				/*
				 * [omega_sr wa omega_ss wb]
				 */
				wt = (omega_sr + wa) / 2.0;
				cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
				if (cosTheta_wt > precision) {
					/*
					 * ANGLES = [omega_sr wa];
					 */
					v_om[0] = omega_sr;
					v_om[1] = wa;
				} else if (cosTheta_wt < -precision) {
					/*
					 * ANGLES = [wa omega_ss];
					 */
					v_om[0] = wa;
					v_om[1] = omega_ss;
				} else {

				}
			} else {
				/*
				 * [wa omega_sr wb omega_ss]
				 */
				wt = (omega_sr + wb) / 2.0;
				cosTheta_wt = cos(wt) * A + sin(wt) * B + C;
				if (cosTheta_wt > precision) {
					/*
					 * ANGLES = [omega_sr wb];
					 */
					v_om[0] = omega_sr;
					v_om[1] = wb;
				} else if (cosTheta_wt < -precision) {
					/*
					 * ANGLES = [wb omega_ss];
					 */
					v_om[0] = wb;
					v_om[1] = omega_ss;
				} else {

				}
			}
		}
	}

	if (v_om[0] == -999.0) {
		*p_nb = 0;
	} else if (v_om[2] == -999.0) {
		*p_nb = 2;
	} else {
		*p_nb = 0;
	}

	/* Patch pour faire fonctionner le code de Mireille */
	v_om[0] = omega_sr;
	v_om[1] = omega_ss;
	*p_nb = 2;

}
