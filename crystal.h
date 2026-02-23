#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <stddef.h>
#include <stdbool.h>
#include "math_helper.h"   

typedef enum { CUBIC, TETRAGONAL, HEXAGONAL, ORTHORHOMBIC, RHOMBOHEDRAL, MONOCLINIC, TRICLINIC } System;
typedef enum { PRIMITIVE, BODY_CENTERED, FACE_CENTERED, BASE_CENTERED} BasisType;


// STRUCTS ------------------------ //

typedef struct {
    Mat3 A;    // conventional cell basis (a,b,c)
    Mat3 B;    // reciprocal basis (a*,b*,c*)
    double a,b,c;
    double alpha,beta,gamma;
    System type;
} Lattice;


typedef struct {
    size_t n;
    Vec3 *pos;   // atomic positions (fractional unit cell) 
    int *Z;
    BasisType type;
} BasisAtoms;


typedef struct {
    HKL hkl;
    double u, v; 
    double intensity; 
} ReciprocalPoint;


typedef struct {
    size_t n;
    ReciprocalPoint *pts; 
    HKL zone; // The normal vector of our plane which slices through the 3D reciprocal space    
} ReciprocalSpace;


typedef struct {
    Lattice lattice;
    BasisAtoms *basis;
    ReciprocalSpace *space;
} Crystal;



void basis_atoms_destroy(BasisAtoms *bas);


bool basis_atoms_resize(size_t n, BasisAtoms *bas);


void rs_destroy(ReciprocalSpace *rs);


bool rs_resize(ReciprocalSpace *rs, size_t n, HKL zone);


void crystal_free(Crystal *crystal);


Crystal* crystal_init(double a, double b, double c, double alpha, double beta, double gamma);


bool rs_basis(Crystal *crystal);


double structure_factor(Crystal *crystal, HKL plane);


bool generate_relp(Crystal *crystal, HKL zone);


bool generate_cell(Crystal *crystal, System sys, BasisType bas);


bool generate_space(Crystal *crystal, System sys, BasisType bas, HKL zone);


bool validate_lat_params(System type, double a, double b, double c, double alpha, double beta, double gamma);


const char *options(System type);


#endif