#include "prop.h"

// calculate phi at x = j
static complex double phij(complex double psijm, complex double psij,
        complex double psijp, double Vj) {
    return psij - I * dt * 0.5 *
        (Vj * psij - (psijm - 2 * psij + psijp) * 0.5 / (dx * dx));
}

void init_x(double* xpos) {
    for (int j = 0; j < nx; j++) {
        xpos[j] = dx * j;
    }
}

void init_psi0(complex double* psi, double* xpos, double wave_num) {
    k = wave_num;
    double sum = 0;
    psi[0] = psi[nx - 1] = 0;
    for (int j = 1; j < nx - 1; j++) {
        psi[j] = exp(-pow(xpos[j] - R * 0.2, 2) * 0.5 / (std * std))
            * cexp(I * k * xpos[j]);
        sum += psi[j] * conj(psi[j]);
    }

    // normalize
    sum = sqrt(sum);
    for (int j = 0; j < nx; j++) {
        psi[j] /= sum;
    }
}

void init_V(double* V, double* xpos) {
    for (int j = 0; j < nx; j++) {
        V[j] = -V0 * exp(-pow((xpos[j] - R * 0.5) / R * a, 26) * 0.5);
    }
}

void init_d(complex double* d, complex double* d2, const double* V) {
    for (int j = 0; j < nx; j++) {
        d2[j] = d[j] = 1 + I * dt * 0.5 * (V[j] + 1 / (dx * dx));
    }
}

void next_phi(const complex double* psi, complex double* phi, const double* V) {
    phi[0] = phij(0, psi[0], psi[1], V[0]);
    for (int j = 1; j < nx - 1; j++) {
        phi[j] = phij(psi[j - 1], psi[j], psi[j + 1], V[j]);
    }
    phi[nx - 1] = phij(psi[nx - 2], psi[nx - 1], 0, V[nx - 1]);
}

// recursive solver
void solver(complex double* psi, complex double* phi, const complex double* d,
        complex double* d2, int j) {
    e = -I * dt / (dx * dx) * 0.25;
    if (j == 1) { // boundary condition (psi(x = 0) = 0)
        solver(psi, phi, d, d2, j + 1);
        psi[j] = (phi[j] - e * psi[j + 1]) / d[j];
    } else if (j == nx - 2) { // boundary condition
        d2[j] = d[j] - e * e / d2[j - 1];
        phi[j] = phi[j] - (e * phi[j - 1] / d2[j - 1]);
        psi[j] = phi[j] / d2[j];
    } else {
        d2[j] = d[j] - e * e / d2[j - 1];
        phi[j] = phi[j] - (e * phi[j - 1] / d2[j - 1]);
        solver(psi, phi, d, d2, j + 1);
        psi[j] = (phi[j] - e * psi[j + 1]) / d2[j];
    }
}

void norm(double* normed, const complex double* psi) {
    for (int j = 0; j < nx; j++) {
        normed[j] = psi[j] * conj(psi[j]);
    }
}
