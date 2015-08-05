/*
 * angle.h
 *
 *  Created on: Apr 16, 2012
 *      Author: Laurent Saboret, Transvalor
 */

#ifndef ANGLE_H_
#define ANGLE_H_

#include <math.h>
#include <float.h>



#ifndef M_PI
#define M_PI 3.141592653589793
#define M_PI_2 1.570796326794897
#endif

#ifdef  __cplusplus
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
} angle_t;

static INLINE angle_t init_angle(double angle) { // angle in radians
    angle_t ths;
    
    ths.angle     = angle;
    
    /* other fields are not yet computed */
    ths.cos_angle = DBL_MAX;
    ths.sin_angle = DBL_MAX;
    ths.tan_angle = DBL_MAX;
    
    return ths;
}

static INLINE angle_t init_angle_cos(double angle, double cos_angle) {
    angle_t ths;
    
    ths.angle     = angle;
    ths.cos_angle = cos_angle;
    
    /* other fields are not yet computed */
    ths.sin_angle = DBL_MAX;
    ths.tan_angle = DBL_MAX;
    
    return ths;
}

static INLINE angle_t init_angle_sin(double angle, double sin_angle) {
    angle_t ths;
    
    ths.angle     = angle;
    ths.sin_angle = sin_angle;
    
    /* other fields are not yet computed */
    ths.cos_angle = DBL_MAX;
    ths.tan_angle = DBL_MAX;
    
    return ths;
}

static INLINE angle_t init_angle_cos_sin(double angle, double cos_angle, double sin_angle) {
    angle_t ths;
    
    ths.angle     = angle;
    ths.cos_angle = cos_angle;
    ths.sin_angle = sin_angle;
    
    /* other fields are not yet computed */
    ths.tan_angle = DBL_MAX;
    
    return ths;
}

// Return -a
static INLINE angle_t minus_angle(const angle_t* a) {
    angle_t ths;
    
    ths.angle     = a->angle;
    ths.cos_angle = a->cos_angle;
    ths.sin_angle = (a->sin_angle != DBL_MAX) ? (-a->sin_angle) : DBL_MAX;
    ths.tan_angle = (a->tan_angle != DBL_MAX) ? (-a->tan_angle) : DBL_MAX;
    
    return ths;
}

// Return PI/2 - a
static INLINE angle_t pi_2_minus_angle(const angle_t* a) {
    angle_t ths;
    
    ths.angle     = M_PI_2 - a->angle;
    ths.cos_angle = a->sin_angle;
    ths.sin_angle = a->cos_angle;
    ths.tan_angle = (a->tan_angle != DBL_MAX) ? (1.0/a->tan_angle) : DBL_MAX;
    
    return ths;
}

// get angle in radians
static INLINE double get_angle(const angle_t* ths) { 
    return ths->angle;
}

// compute cosinus only once
static INLINE double get_cos(const angle_t* ths) { 
    if (ths->cos_angle == DBL_MAX)
        ((angle_t*)ths)->cos_angle = cos(ths->angle);
    
    return ths->cos_angle;
}

// compute sinus only once
static INLINE double get_sin(const angle_t* ths) { 
    if (ths->sin_angle == DBL_MAX)
        ((angle_t*)ths)->sin_angle = sin(ths->angle);
    
    return ths->sin_angle;
}

// compute tangent only once
static INLINE double get_tan(const angle_t* ths) { 
    if (ths->tan_angle == DBL_MAX) {
        if (ths->cos_angle != DBL_MAX && ths->sin_angle != DBL_MAX)
            ((angle_t*)ths)->tan_angle = ths->sin_angle / ths->cos_angle;
        else
            ((angle_t*)ths)->tan_angle = tan(ths->angle);
    }
    
    return ths->tan_angle;
}

// force cosinus, sinus and tangent computation
static INLINE void fill_cache(const angle_t* ths) { 
    (void) get_cos(ths);
    (void) get_sin(ths);
    (void) get_tan(ths);
}

#ifdef __cplusplus
}
#endif

#endif /* ANGLE_H_ */
