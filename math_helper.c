#include "math_helper.h"


// METHODS ------------------------ //

// Get a normal Vec3 to input Vec3 n
Vec3 v3_unit_normal(Vec3 n) {
    Vec3 u;
    if (fabs(n.x) <= fabs(n.y) && fabs(n.x) <= fabs(n.z)) { 
		u = (Vec3){1,0,0}; 
	}
    else if (fabs(n.y) <= fabs(n.z)) { 
		u = (Vec3){0,1,0}; 
	}
    else { 
		u = (Vec3){0,0,1}; 
	}
    return v3_normalize(v3_cross(n, u));
}


// // Use Gram-Schmidt to orthonormalize a basis, represented by Mat3 A
// bool gram_schmidt(const Mat3 A, Mat3 *Out) {
// 	const double eps  = 1e-12;
// 	const double eps2 = eps * eps;
	
// 	Vec3 v1 = mat3_col(A, 0);
// 	Vec3 v2 = mat3_col(A, 1);
// 	Vec3 v3 = mat3_col(A, 2);

// 	// e1 = v1 / ||v1||
// 	double n1 = v3_dot(v1, v1);
// 	if (n1 < eps2) { return false; }
// 	Vec3 e1 = v3_scale(v1, 1.0 / sqrt(n1));

// 	// u2 = v2 - proj of v2 on e1
// 	Vec3 u2 = v3_sub(v2, v3_proj(v2, e1));
// 	double n2 = v3_dot(u2, u2);
// 	if (n2 < eps2) { return false; }
// 	Vec3 e2 = v3_scale(u2, 1.0 / sqrt(n2));

// 	// u3 = v3 with components along e1 and e2 removed 
// 	Vec3 u3 = v3_sub(v3, v3_proj(v3, e1));
// 	u3 = v3_sub(u3, v3_proj(u3, e2));
// 	double n3 = v3_dot(u3, u3);
// 	if (n3 < eps2) { return false; }
// 	Vec3 e3 = v3_scale(u3, 1.0 / sqrt(n3));

// 	*Out = v3_to_mat3(e1, e2, e3);
// 	return true;
// };


