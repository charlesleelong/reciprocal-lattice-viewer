#ifndef MATH_HELPER_H
#define MATH_HELPER_H

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef DEG2RAD
#define DEG2RAD (PI / 180.0)
#endif


// Contains structs/methods related to useful linear algebra operations

// STRUCTS ------------------------ //

typedef struct {
	double x, y, z;
} Vec3;


typedef struct {
	double M[3][3]; 
} Mat3;


typedef struct { 
	int h,k,l; 
} HKL;


// METHODS ------------------------ //

// Construct new Vec3
static inline Vec3 v3(double x, double y, double z);


// Multiply Vec3 by scalar a
static inline Vec3 v3_scale(const Vec3 vec, double a) {
	return (Vec3){
		vec.x * a,
		vec.y * a,
		vec.z * a
	};
};


// Add two Vec3
static inline Vec3 v3_add(const Vec3 a, const Vec3 b) {
	return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}


// Subtract two Vec3
static inline Vec3 v3_sub(const Vec3 a, const Vec3 b) {
	return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}


// Dot product of two Vec3
static inline double v3_dot(const Vec3 a, const Vec3 b) {
	return ( a.x * b.x ) + ( a.y * b.y ) + ( a.z * b.z );
};


// Cross product of two Vec3
static inline Vec3 v3_cross(const Vec3 a, const Vec3 b) {
	return (Vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
};


static inline double v3_magnitude(Vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}


static inline Vec3 v3_normalize(Vec3 v) {
    double n2 = v3_dot(v,v);
    if (n2 < 1e-24) return (Vec3){0,0,0};
    return v3_scale(v, 1.0/sqrt(n2));
};


static inline Vec3 v3_project(Vec3 v, Vec3 u) {
    double uu = v3_dot(u, u);
    if (uu < 1e-24) return (Vec3){0,0,0};
    return v3_scale(u, v3_dot(v,u) / uu);
}


// Get column vector "col" from Mat3 A
static inline Vec3 mat3_col(Mat3 A, int col) {
    return (Vec3){
        A.M[0][col],
        A.M[1][col],
        A.M[2][col]
    };
};


// Get row vect "row" from Mat3 A
static inline Vec3 mat3_row(Mat3 A, int row) {
    return (Vec3){
        A.M[row][0],
        A.M[row][1],
        A.M[row][2]
    };
};


// Convert 3 Vec3 into Mat3 (column vectors -> matrix)
static inline Mat3 v3_to_mat3(Vec3 a, Vec3 b, Vec3 c) {
	return (Mat3){
        .M = {
            { a.x, b.x, c.x },
            { a.y, b.y, c.y },
            { a.z, b.z, c.z } 
        }
    };
};


Vec3 v3_unit_normal(Vec3 n);


static inline HKL hkl_scale(HKL a, int s)
{
    return (HKL){ a.h*s, a.k*s, a.l*s };
}


static inline HKL hkl_add(HKL a, HKL b)
{
    return (HKL){ a.h + b.h, a.k + b.k, a.l + b.l };
}


static inline double hkl_dot(HKL a, HKL b)
{
    return a.h * b.h + a.k * b.k + a.l * b.l;
}


static inline HKL hkl_cross(HKL a, HKL b)
{
    return (HKL){
        a.k*b.l - a.l*b.k,
        a.l*b.h - a.h*b.l,
        a.h*b.k - a.k*b.h
    };
}


static inline HKL hkl_normal(HKL a) {
    if (fabs(a.h) < fabs(a.k) && fabs(a.h) < fabs(a.l)) {
        return (HKL) {1, 0, 0};
    }
    else if (fabs(a.k) < fabs(a.h) && fabs(a.k) < fabs(a.l)) {
        return (HKL) {0, 1, 0};
    }
    else {
        return (HKL) {0, 0, 1};
    }
}


static inline Vec3 hkl_to_v3(HKL a) {
    return (Vec3) {(double)a.h, (double)a.k, (double)a.l};
}


#endif