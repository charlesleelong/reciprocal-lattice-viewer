#include "crystal.h"

#include <stdlib.h>
#include <math.h>
#include <complex.h>


void basis_atoms_destroy(BasisAtoms *bas) {
    if (!bas) { return; }

    free(bas->pos);
    free(bas->Z);
    free(bas);

    return;
}


bool basis_atoms_resize(size_t n, BasisAtoms *bas) {
    if (!bas) { return false; }

    if (bas->n == n) { return true; }

    if (n == 0) {
        free(bas->pos);
        free(bas->Z);
        bas->pos = NULL;
        bas->Z = NULL;
        bas->n = 0;
        return true; 
    }

    Vec3 *new_pos = malloc(n * sizeof(*new_pos));
    int *new_Z = malloc(n * sizeof(*new_Z));

    if (!new_pos || !new_Z) {
        free(new_pos);
        free(new_Z);
        return false;
    }

    free(bas->pos);
    free(bas->Z);
    bas->pos = new_pos;
    bas->Z = new_Z;
    bas->n = n;

    return true;
}


void rs_destroy(ReciprocalSpace *rs) {
    if (!rs) { return; }

    free(rs->pts);
    free(rs);

    return;
}


bool rs_resize(ReciprocalSpace *rs, size_t n, HKL zone) {
    if (!rs) { return false; }

    rs->zone = zone;

    if (rs->n == n) {
        return true;
    }

    if (n == 0) {
        free(rs->pts);
        rs->pts = NULL;
        rs->n = 0;
        return true;
    }

    ReciprocalPoint *new_pts = calloc(n, sizeof(*new_pts));
    if (!new_pts) {
        return false; 
    }

    free(rs->pts);
    rs->pts = new_pts;
    rs->n = n;

    return true;
}


void crystal_free(Crystal *crystal) {
    if (!crystal) { return; }

    basis_atoms_destroy(crystal->basis);
    rs_destroy(crystal->space); 
    free(crystal);

    return;
}


Crystal* crystal_init(double a, double b, double c, double alpha, double beta, double gamma) {
    Crystal *crystal = calloc(1, sizeof(*crystal));
    if (!crystal) { return NULL; }

    crystal->basis = calloc(1, sizeof(*crystal->basis));
    if (!crystal->basis) { crystal_free(crystal); return NULL; }

    crystal->space = calloc(1, sizeof(*crystal->space));
    if (!crystal->space) { crystal_free(crystal); return NULL; }

    crystal->lattice.a = a;
    crystal->lattice.b = b;
    crystal->lattice.c = c;
    crystal->lattice.alpha = alpha;
    crystal->lattice.beta = beta;
    crystal->lattice.gamma = gamma;

    return crystal;
}


// METHODS ------------------------ //

// Construct reciprocal vector basis from conventional basis
bool rs_basis(Crystal *crystal) {
	Vec3 a = mat3_col(crystal->lattice.A, 0);
	Vec3 b = mat3_col(crystal->lattice.A, 1);
	Vec3 c = mat3_col(crystal->lattice.A, 2);

	double vol = v3_dot(a, v3_cross(b, c)); // Cell volume
	double scalar = (2 * PI) / vol; // Physics convention introduces factor of 2PI to numerator

	Vec3 a_r = v3_scale(v3_cross(b, c), scalar);
	Vec3 b_r = v3_scale(v3_cross(c, a), scalar);
	Vec3 c_r  = v3_scale(v3_cross(a, b), scalar);

	crystal->lattice.B = v3_to_mat3(a_r, b_r, c_r);
	return true;	
};


double structure_factor(Crystal *crystal, HKL plane) {
    double complex F = 0 + 0 * I; 
    Vec3 atom;
    int Z;
    double f = 1.0;
    for(size_t i = 0; i < crystal->basis->n; i++) {
        atom = crystal->basis->pos[i];
        double phase = 2 * PI * v3_dot(atom, hkl_to_v3(plane));
        // if (crystal->basis->Z) {
        //     f = (double)crystal->basis->Z[i];    // scattering factor (approximated as atomic number)
        // }
        F += f * cexp(I * phase); 
    }  

    return pow(cabs(F), 2);
}


// Generate reciprocal lattice points from a Crystal struct chosen plane normal
bool generate_relp(Crystal *crystal, HKL zone) {
    // Normal vector cannot be zero
    if ( !crystal || (zone.h == 0 && zone.k == 0 && zone.l == 0) ) {
        return false; 
    }

    // Reciprocal vectors 
    Vec3 b1 = mat3_col(crystal->lattice.B, 0);    
    Vec3 b2 = mat3_col(crystal->lattice.B, 1);    
    Vec3 b3 = mat3_col(crystal->lattice.B, 2);    

    Vec3 zone_rs = v3_add(
        v3_add(v3_scale(b1, (double)zone.h),
               v3_scale(b2, (double)zone.k)),
               v3_scale(b3, (double)zone.l)
    );

    zone_rs = v3_normalize(zone_rs);

    int z_h = zone.h;
    int z_k = zone.k;
    int z_l = zone.l;

    Vec3 q;

    int H = 15;
    int h,k,l;
    int count = 0;
    for (h = -H; h <= H; h++) {
        for (k = -H; k <= H; k++) {
            for (l = -H; l <= H; l++) {
                if (h * z_h + k * z_k + l * z_l == 0) {
                    count++;
                }
            }
        }
    }

    if (!rs_resize(crystal->space, count, zone)) { return false; }

    Vec3 e1 = v3_unit_normal(zone_rs);
    Vec3 e2 = v3_cross(zone_rs, e1);

    Vec3 x_dir = (Vec3){1, 0, 0};
    if (v3_dot(e1, x_dir) < 0) {
        e1 = v3_scale(e1, -1.0);
        e2 = v3_scale(e2, -1.0);  
    }

    int n = 0;
    for (h = -H; h <= H; h++) {
        for (k = -H; k <= H; k++) {
            for (l = -H; l <= H; l++) {
                if (h * z_h + k * z_k + l * z_l == 0) {
                    q = v3_add(
                        v3_add(v3_scale(b1, (double)h),
                           v3_scale(b2, (double)k)),
                           v3_scale(b3, (double)l)
                    );
                    HKL plane = (HKL){h, k, l};
                    crystal->space->pts[n].hkl = plane; 
                    crystal->space->pts[n].u = v3_dot(q, e1);
                    crystal->space->pts[n].v = v3_dot(q, e2);
                    crystal->space->pts[n].intensity = structure_factor(crystal, plane);
                    n++;
                }
            }
        }
    }
    
    return true;
}


bool generate_cell(Crystal *crystal, System sys, BasisType bas) {
    if (!crystal || !crystal->basis) return false;

    double a = crystal->lattice.a;
    double b = crystal->lattice.b;
    double c = crystal->lattice.c;
    double alpha = crystal->lattice.alpha;
    double beta = crystal->lattice.beta;
    double gamma = crystal->lattice.gamma;

    switch(sys) {
        case CUBIC: 
            switch(bas) {
                case PRIMITIVE:
                case BODY_CENTERED:
                case FACE_CENTERED:
                    crystal->lattice.A = (Mat3){ 
                        .M = {
                                { a,    0,    0 },
                                { 0,    a,    0 },
                                { 0,    0,    a }
                    }};
                    break;
                default: return false;
            }
            break;
        case TETRAGONAL:
            switch(bas) {
                case PRIMITIVE:
                case BODY_CENTERED:
                    crystal->lattice.A = (Mat3){ 
                        .M = {
                                { a,    0,    0 },
                                { 0,    a,    0 },
                                { 0,    0,    c }
                    }};
                    break;
                default: return false;
            }
            break;
        case HEXAGONAL:
            switch(bas) {
                case PRIMITIVE:
                    crystal->lattice.A = (Mat3){ 
                        .M = {
                            { a,    -a / 2,           0 },
                            { 0,    sqrt(3) * a / 2,  0 },  
                            { 0,    0,                c }
                        }
                    };
                    break;
                default: return false;
            }
            break;
        case ORTHORHOMBIC:
            switch(bas) {
                case PRIMITIVE:
                case BODY_CENTERED:
                case FACE_CENTERED:
                case BASE_CENTERED:
                    crystal->lattice.A = (Mat3){ 
                        .M = {
                                { a,    0,    0 },
                                { 0,    b,    0 },
                                { 0,    0,    c }
                    }};
                    break;
                default: return false;
            }
            break;
        case MONOCLINIC:
            switch(bas) {
                case PRIMITIVE:
                case BASE_CENTERED:
                    crystal->lattice.A = (Mat3){ 
                        .M = {
                                { a,    0,    c * cos((DEG2RAD * beta)) },
                                { 0,    b,                            0 },
                                { 0,    0,    c * sin((DEG2RAD * beta)) }
                    }};
                    break;
                default: return false;
            }
            break;
        case TRICLINIC:
        case RHOMBOHEDRAL: {
            switch(bas) {
                case PRIMITIVE:
                    double c_x = c * cos((DEG2RAD * beta));
                    double c_y = c * (cos((DEG2RAD * alpha)) - cos((DEG2RAD * beta)) * cos((DEG2RAD * gamma))) / sin((DEG2RAD * gamma));
                    double c_z = sqrt(c * c - c_x * c_x - c_y * c_y);
                    crystal->lattice.A = (Mat3){ 
                        .M = {
                                { a,    b * cos((DEG2RAD * gamma)),    c_x},
                                { 0,    b * sin((DEG2RAD * gamma)),    c_y},
                                { 0,    0,                             c_z}
                    }};
                    break;
                default: return false;
            }
            break;
        }
        default: return false;
    }
    switch(bas) {
        case PRIMITIVE:
            if (sys == RHOMBOHEDRAL) {
                if (!basis_atoms_resize(3, crystal->basis)) { return false; }
                crystal->basis->pos[0] = (Vec3){0.0, 0.0, 0.0};
                crystal->basis->pos[1] = (Vec3){2.0/3.0, 1.0/3.0, 1.0/3.0};
                crystal->basis->pos[2] = (Vec3){1.0/3.0, 2.0/3.0, 2.0/3.0};
                break;
            }
            if (!basis_atoms_resize(1, crystal->basis)) { return false; }
            crystal->basis->pos[0] = (Vec3){0, 0, 0};
            break;
        case BODY_CENTERED:
            if (!basis_atoms_resize(2, crystal->basis)) { return false; }
            crystal->basis->pos[0] = (Vec3){0, 0, 0};
            crystal->basis->pos[1] = (Vec3){0.5, 0.5, 0.5};
            break;
        case FACE_CENTERED:
            if (!basis_atoms_resize(4, crystal->basis)) { return false; }
            crystal->basis->pos[0] = (Vec3){0, 0, 0};
            crystal->basis->pos[1] = (Vec3){0.5, 0.5, 0};
            crystal->basis->pos[2] = (Vec3){0.5, 0, 0.5};
            crystal->basis->pos[3] = (Vec3){0, 0.5, 0.5};
            break;
        case BASE_CENTERED:
            if (!basis_atoms_resize(2, crystal->basis)) { return false; }
            crystal->basis->pos[0] = (Vec3){0, 0, 0};
            crystal->basis->pos[1] = (Vec3){0.5, 0.5, 0};
            break;
        default: return false;
    }

    crystal->lattice.a = a; crystal->lattice.b = b; crystal->lattice.c = c;
    crystal->lattice.alpha = alpha; crystal->lattice.beta = beta; crystal->lattice.gamma = gamma;
    crystal->lattice.type = sys;
    crystal->basis->type = bas;

    return true;
}


bool generate_space(Crystal *crystal, System sys, BasisType bas, HKL zone) {
    if (!crystal || !crystal->basis || !crystal->space) return false;

    crystal->lattice.type = sys;
    crystal->basis->type = bas;

    if (!generate_cell(crystal, sys, bas)) return false;
    if (!rs_basis(crystal)) return false;
    if (!generate_relp(crystal, zone)) return false;
    
    crystal->space->zone = zone;

    return true;
}


bool validate_lat_params(System type, double a, double b, double c, double alpha, double beta, double gamma) {
    switch (type) {
        case CUBIC:
            return (a == b) && (b == c) &&
                   (alpha == 90.0) && (beta == 90.0) && (gamma == 90.0);
        case TETRAGONAL:
            return (a == b) &&
                   (alpha == 90.0) && (beta == 90.0) && (gamma == 90.0);
        case HEXAGONAL:
            return (a == b) &&
                   (alpha == 90.0) && (beta == 90.0) && (gamma == 120.0);
        case ORTHORHOMBIC:
            return (alpha == 90.0) && (beta == 90.0) && (gamma == 90.0);
        case RHOMBOHEDRAL: 
            return (a == b) && (b == c) &&
                   (alpha == beta) && (beta == gamma);
        case MONOCLINIC:
            return (alpha == 90.0) && (gamma == 90.0);
        case TRICLINIC:
            return true;
        default:
            return false;
    }
}


const char *options(System type) {
    switch (type) {
        case CUBIC:
            return "PRIMITIVE;BODY_CENTERED;FACE_CENTERED";
            break;
        case TETRAGONAL:
            return "PRIMITIVE;BODY_CENTERED";
            break;
        case HEXAGONAL:
            return "PRIMITIVE";
            break;
        case ORTHORHOMBIC:
            return "PRIMITIVE;BODY_CENTERED;FACE_CENTERED;BASE_CENTERED";
            break;
        case RHOMBOHEDRAL:  
            return "PRIMITIVE";
            break;
        case MONOCLINIC:
            return "PRIMITIVE;BASE_CENTERED";
            break;
        case TRICLINIC:
            return "PRIMITIVE";
            break;
        default: return "INVALID CRYSTAL SYSTEM";
    }
} 