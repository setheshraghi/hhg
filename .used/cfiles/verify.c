#include <stdlib.h>
#include "prop.h"

// calculates the existential probabiltiy of the particle
static double exist_prob(const complex double* psi1,
        const complex double* psi2) {
    double sum = 0;
    for (int j = 0; j < nx; j++) {
        sum += conj(psi1[j]) * psi2[j];
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
    complex double* psi0 = (complex double*) calloc(nx, sizeof(complex double));
    double* normed = (double*) calloc(nx, sizeof(double));
    complex double* phi = (complex double*) calloc(nx, sizeof(complex double));
    double* V = (double*) calloc(nx, sizeof(double));
    // diagonal vectors, d2 is used by the solver function
    complex double* d = (complex double*) calloc(nx, sizeof(complex double));
    complex double* d2 = (complex double*) calloc(nx, sizeof(complex double));
    complex double* psi0_tilda = (complex double*)
        calloc(nx, sizeof(complex double));

    float k, dt;
    scanf("k=%f,dt=%f", &k, &dt);
    init_x(xpos);
    init_psi0(psi0, xpos, k);
    norm(normed, psi0);
    fprintf(results, "%lf %lf %lf %lf %lf\n", 0.0, exist_prob(psi0, psi0),
            expect_x(normed), expect_H(psi0, V), expect_p(psi0));

    init_V(V, xpos);
    init_d(d, d2, V, dt);
    for (int j = 0; j < nx; j++) {
        psi0_tilda[j] = psi0[j];
    }

    for (int n = 1; n < nt; n++) {
        next_phi(psi0_tilda, phi, V, dt);
        solver(psi0_tilda, phi, d, d2, 1, dt);
        norm(normed, psi0_tilda);
        fprintf(results, "%.2lf %lf %lf %lf %lf\n",
                n * dt,
                exist_prob(psi0_tilda, psi0_tilda),
                expect_x(normed),
                expect_H(psi0_tilda, V),
                expect_p(psi0_tilda));
    }

    dt = -dt;
    init_d(d, d2, V, dt);
    for (int n = 1; n < nt; n++) {
        next_phi(psi0_tilda, phi, V, dt);
        solver(psi0_tilda, phi, d, d2, 1, dt);
    }
    printf("\t autocorr norm: %lf\n", exist_prob(psi0, psi0_tilda));

    free(xpos);
    free(psi0);
    free(psi0_tilda);
    free(phi);
    free(V);
    free(normed);
    free(d);
    free(d2);
    fclose(results);
}
