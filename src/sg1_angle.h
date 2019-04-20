/*
 * angle.h
 *
 *  Created on: Apr 16, 2012
 *      Author: Laurent Saboret, Transvalor
 */

#ifndef SG1_ANGLE_H_
#define SG1_ANGLE_H_

#ifdef __cplusplus

#include <cmath>
#include <cassert>

// same as isnan()
inline static bool _sg1_isnan(double x) { return std::isnan(x); }

#else

#include <math.h>
#include <assert.h>

inline static int _sg1_isnan(double x) { return isnan(x); }

#endif

/* M_PI is not part of C99 or C++11 standard, thus define ours. */
#ifdef M_PI
#define _SG1_PI     M_PI
#define _SG1_PI_2   M_PI_2
#else
#define _SG1_PI     3.141592653589793
#define _SG1_PI_2   1.570796326794897
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/* This structure represents an angle expressed in radians.
 * Calls to cos(), sin(), tan() are cached.
 * 
 * WARNING: always use accessors to read/write this structure's fields!
 */
typedef struct {
/*private: */
                double angle; // angle in radians
    /*mutable*/ double cos_angle; // cos(angle)
    /*mutable*/ double sin_angle; // sin(angle)
    /*mutable*/ double tan_angle; // tan(angle)
} sg1_angle;

inline static sg1_angle sg1_init_angle(double angle) { // angle in radians
    sg1_angle ths;
    
    assert(!_sg1_isnan(angle));
    
    ths.angle     = angle;
    
    /* other fields are not yet computed */
    ths.cos_angle = NAN;
    ths.sin_angle = NAN;
    ths.tan_angle = NAN;
    
    return ths;
}

inline static sg1_angle sg1_init_angle_cos(double angle, double cos_angle) {
    sg1_angle ths;
    
    assert(!_sg1_isnan(angle));
    
    ths.angle     = angle;
    ths.cos_angle = cos_angle;
    
    /* other fields are not yet computed */
    ths.sin_angle = NAN;
    ths.tan_angle = NAN;
    
    return ths;
}

inline static sg1_angle sg1_init_angle_sin(double angle, double sin_angle) {
    sg1_angle ths;
    
    assert(!_sg1_isnan(angle));
    
    ths.angle     = angle;
    ths.sin_angle = sin_angle;
    
    /* other fields are not yet computed */
    ths.cos_angle = NAN;
    ths.tan_angle = NAN;
    
    return ths;
}

inline static sg1_angle sg1_init_angle_tan(double angle, double tan_angle) {
    sg1_angle ths;
    
    assert(!_sg1_isnan(angle));
    
    ths.angle     = angle;
    ths.tan_angle = tan_angle;
    
    /* other fields are not yet computed */
    ths.cos_angle = NAN;
    ths.sin_angle = NAN;
    
    return ths;
}

inline static sg1_angle sg1_init_angle_cos_sin(double angle, double cos_angle, double sin_angle) {
    sg1_angle ths;
    
    assert(!_sg1_isnan(angle));
    
    ths.angle     = angle;
    ths.cos_angle = cos_angle;
    ths.sin_angle = sin_angle;
    
    /* other fields are not yet computed */
    ths.tan_angle = NAN;
    
    return ths;
}

// Return -a
inline static sg1_angle sg1_minus_angle(const sg1_angle* a) {
    sg1_angle ths;
    
    ths.angle     = -a->angle;
    ths.cos_angle = a->cos_angle;
    ths.sin_angle = -a->sin_angle;
    ths.tan_angle = -a->tan_angle;
    
    return ths;
}

// Return PI/2 - a
inline static sg1_angle sg1_pi_2_minus_angle(const sg1_angle* a) {
    sg1_angle ths;
    
    ths.angle     = _SG1_PI_2 - a->angle;
    ths.cos_angle = a->sin_angle;
    ths.sin_angle = a->cos_angle;
    ths.tan_angle = 1.0/a->tan_angle;
    
    return ths;
}

// Return a + PI 
inline static sg1_angle sg1_angle_plus_pi(const sg1_angle* a) {
    sg1_angle ths;
    
    ths.angle     = a->angle + _SG1_PI;
    ths.cos_angle = -a->cos_angle;
    ths.sin_angle = -a->sin_angle;
    ths.tan_angle = a->tan_angle;
    
    return ths;
}

// Return a + 2*PI 
inline static sg1_angle sg1_angle_plus_2pi(const sg1_angle* a) {
    sg1_angle ths;
    
    ths.angle     = a->angle + 2.0*_SG1_PI;
    ths.cos_angle = a->cos_angle;
    ths.sin_angle = a->sin_angle;
    ths.tan_angle = a->tan_angle;
    
    return ths;
}

// Return a - 2*PI 
inline static sg1_angle sg1_angle_minus_2pi(const sg1_angle* a) {
    sg1_angle ths;
    
    ths.angle     = a->angle - 2.0*_SG1_PI;
    ths.cos_angle = a->cos_angle;
    ths.sin_angle = a->sin_angle;
    ths.tan_angle = a->tan_angle;
    
    return ths;
}

// get angle in radians
inline static double sg1_get_angle(const sg1_angle* ths) {
    return ths->angle;
}

// compute cosinus (cached)
inline static double sg1_get_cos(const sg1_angle* ths) {
    if (_sg1_isnan(ths->cos_angle))
        ((sg1_angle*)ths)->cos_angle = cos(ths->angle);
    
    return ths->cos_angle;
}

// compute sinus (cached)
inline static double sg1_get_sin(const sg1_angle* ths) {
    if (_sg1_isnan(ths->sin_angle))
        ((sg1_angle*)ths)->sin_angle = sin(ths->angle);
    
    return ths->sin_angle;
}

// compute tangent (cached)
inline static double sg1_get_tan(const sg1_angle* ths) {
    if (_sg1_isnan(ths->tan_angle))
        ((sg1_angle*)ths)->tan_angle = tan(ths->angle);
    
    return ths->tan_angle;
}

// force cosinus, sinus and tangent computation
inline static void sg1_fill_cache(const sg1_angle* ths) {
    (void) sg1_get_cos(ths);
    (void) sg1_get_sin(ths);
    (void) sg1_get_tan(ths);
}


#ifdef __cplusplus
}
#endif

#endif /* SG1_ANGLE_H_ */
