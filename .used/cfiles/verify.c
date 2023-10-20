#include <stdlib.h>
#include "prop.h"

// calculates the existential probabiltiy of the particle
static double exist_prob(const complex double* psi1,
        const complex double* psi2) {
    double sum = 0;
    for (int j = 0; j < nx; j++) {
        sum += psi1[j] * conj(psi2[j]);
    }
    return sum;
}

// Hamiltonian operator at x = j
static complex double hamiltonian(complex double jm, complex double j,
        complex double jp, double Vj) {
    return Vj * j - (jm - 2 * j + jp) * 0.5 / (dx * dx);
}

// momentum operator at x = jm + 1
static complex double momentum(complex double jm, complex double jp) {
    return -I * (jp - jm) * 0.5 / dx;
}

// return expectation value of x position
static double expect_x(const double* normed) {
    double sum = 0;
    for (int j = 0; j < nx; j++) {
        sum += normed[j] * j * dx;
    }
    return sum;
}

// return expectation value of Hamiltonian
static double expect_H(const complex double* psi, const double* V) {
    double sum = 0;
    sum += conj(psi[0]) * hamiltonian(0, psi[0], psi[1], V[0]);
    for (int j = 1; j < nx - 1; j++) {
        sum += conj(psi[j]) * hamiltonian(psi[j - 1], psi[j], psi[j + 1], V[j]);
    }
    sum += conj(psi[nx - 1]) *
        hamiltonian(psi[nx - 2], psi[nx - 1], 0, V[nx - 1]);
    return sum;
}

// return expectation value of momentum
static double expect_p(const complex double* psi) {
    double sum = 0;
    sum += conj(psi[0]) * momentum(0, psi[1]);
    for (int j = 1; j < nx - 1; j++) {
        sum += conj(psi[j]) * momentum(psi[j - 1], psi[j + 1]);
    }
    sum += conj(psi[nx - 1]) * momentum(psi[nx - 2], 0);
    return sum;
}

int main() {
    FILE* results = fopen(".used/datafiles/verify.txt", "w");
    double* xpos = (double*) calloc(nx, sizeof(double));
    complex double* psi = (complex double*) calloc(nx, sizeof(complex double));
    double* normed = (double*) calloc(nx, sizeof(double));
    complex double* phi = (complex double*) calloc(nx, sizeof(complex double));
    double* V = (double*) calloc(nx, sizeof(double));
    // diagonal vectors, d2 is used by the solver function
    complex double* d = (complex double*) calloc(nx, sizeof(complex double));
    complex double* d2 = (complex double*) calloc(nx, sizeof(complex double));
    complex double* psi_tilda = (complex double*)
        calloc(nx, sizeof(complex double));

    init_x(xpos);
    init_psi0(psi, xpos, 3);
    norm(normed, psi);
    fprintf(results, "%lf %lf %lf %lf %lf\n", 0.0, exist_prob(psi, psi),
            expect_x(normed), expect_H(psi, V), expect_p(psi));

    init_V(V, xpos);
    init_d(d, d2, V);

    for (int n = 1; n < nt; n++) {
        next_phi(psi, phi, V);
        solver(psi, phi, d, d2, 1);
        norm(normed, psi);
        fprintf(results, "%.2lf %lf %lf %lf %lf\n",
                n * dt,
                exist_prob(psi, psi),
                expect_x(normed),
                expect_H(psi, V),
                expect_p(psi));
    }

    free(xpos);
    free(psi);
    free(psi_tilda);
    free(phi);
    free(V);
    free(normed);
    free(d);
    free(d2);
    fclose(results);
}
