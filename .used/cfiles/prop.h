#ifndef PROP_H
#define PROP_H
#include <stdio.h>
#include <math.h>
#include <complex.h>

static double k; // initialized in init_psi0
static const int nt = 2222;
static const double dx = 0.2;
static const int nx = 2000;
static const double R = dx * nx;
static const double std = 10;
static const double V0 = 1;
static const double a = R * 0.005;
static complex double e; // value is reassigned in the solver function

// initializes all x values
void init_x(double* xpos);

// initializes psi at t = 0
void init_psi0(complex double* psi, double* xpos, double wave_num);

// initializes potential energy function
void init_V(double* V, double* xpos);

// initializes the diagonal of the matrix used to propogate
void init_d(complex double* d, complex double* d2, const double* V, double dt);

// calculate the next phi
void next_phi(const complex double* psi, complex double* phi, const double* V,
double dt);

// solves linear system where the coefficient matrix is tridiagonal
void solver(complex double* psi, complex double* phi, const complex double* d,
        complex double* d2, int j, double dt);

// finds norm of psi
void norm(double* normed, const complex double* psi);

#endif
