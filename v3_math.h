#ifndef "math_helper.h"
#define "math_helper.h"


// STRUCTS ------------------------ //

typedef struct {
	double x, y, z;
} Vec3;


typedef struct {
	double M[3][3]; 
} Mat3;


typedef struct { int h,k,l; } HKL;


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


static inline Vec3 v3_norm(Vec3 v) {
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

// REMOVE
// Project Vec3 v onto Vec3 u
static inline Vec3 v3_proj(Vec3 v, Vec3 u) {
	double uu = v3_dot(u, u);
	if (fabs(uu) < 1e-12) {  // Assume u is 0-vector in this case, so return it
		return (Vec3){0, 0, 0}; 
	}
	return v3_scale(u, v3_dot(v, u) / uu);
};

// REMOVE
// Get a normal Vec3 to input Vec3 n
static inline Vec3 v3_normal(Vec3 n) {
	Vec3 u;
	// Finds "least-parallel" axis vector
	if (fabs(n.x) < fabs(n.y) && fabs(n.x) < fabs(n.z)) {
        u = (Vec3){1, 0, 0};
    } else if (fabs(n.y) < fabs(n.z)) {
        u = (Vec3){0, 1, 0};
    } else {
        u = (Vec3){0, 0, 1};
    }
    return v3_cross(n, u);
}


// Use Gram-Schmidt to orthonormalize a basis, represented by Mat3 A
static inline bool gram_schmidt(const Mat3 A, Mat3 *Out) {
	const double eps  = 1e-12;
	const double eps2 = eps * eps;
	
	Vec3 v1 = mat3_col(A, 0);
	Vec3 v2 = mat3_col(A, 1);
	Vec3 v3 = mat3_col(A, 2);

	// e1 = v1 / ||v1||
	double n1 = v3_dot(v1, v1);
	if (n1 < eps2) { return false; }
	Vec3 e1 = v3_scale(v1, 1.0 / sqrt(n1));

	// u2 = v2 - proj of v2 on e1
	Vec3 u2 = v3_sub(v2, v3_proj(v2, e1));
	double n2 = v3_dot(u2, u2);
	if (n2 < eps2) { return false; }
	Vec3 e2 = v3_scale(u2, 1.0 / sqrt(n2));

	// u3 = v3 with components along e1 and e2 removed 
	Vec3 u3 = v3_sub(v3, v3_proj(v3, e1));
	u3 = v3_sub(u3, v3_proj(u3, e2));
	double n3 = v3_dot(u3, u3);
	if (n3 < eps2) { return false; }
	Vec3 e3 = v3_scale(u3, 1.0 / sqrt(n3));

	*Out = v3_to_mat3(e1, e2, e3);
	return true;
};


static inline HKL hkl_scale(HKL a, int s)
{
    return (HKL){ a.h*s, a.k*s, a.l*s };
}


static inline HKL hkl_cross(HKL a, HKL b)
{
    return (HKL){
        a.k*b.l - a.l*b.k,
        a.l*b.h - a.h*b.l,
        a.h*b.k - a.k*b.h
    };
}


static inline HKL hkl_add(HKL a, HKL b)
{
    return (HKL){ a.h + b.h, a.k + b.k, a.l + b.l };
}


#endif