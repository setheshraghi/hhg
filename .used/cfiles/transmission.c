#include <stdlib.h>
#include "prop.h"

static double T(double k) {
    double E = k * k * 0.5;
    return 1 / (1 + (V0 * V0 * 0.25 / (E * (E + V0)))
            * pow(sin(2 * a * sqrt(2 * (E + V0))), 2));
}

static void Fprint(FILE* data, complex double* psi, double n, double error) {
    double F = 0;
    for (int j = nx * 0.55; j < nx; j++) {
        F += psi[j] * conj(psi[j]);
    }
    fprintf(data, "%lf %lf %lf %lf\n", n, F, error, T(n)); // A = 1 => T = F
}

static void Tprint(FILE* data, complex double* psi, double n) {
    double F = 0;
    for (int j = nx * 0.55; j < nx; j++) {
        F += psi[j] * conj(psi[j]);
    }
    fprintf(data, "%lf %lf\n", n, F);
}

int main() {
    FILE* trans_data = fopen(".used/datafiles/transmission.dat", "w");
    FILE* F_data = fopen(".used/datafiles/F.dat", "w");
    double* xpos = (double*) calloc(nx, sizeof(double));
    complex double* psi = (complex double*) calloc(nx, sizeof(complex double));
    double* normed = (double*) calloc(nx, sizeof(double));
    complex double* phi = (complex double*) calloc(nx, sizeof(complex double));
    double* V = (double*) calloc(nx, sizeof(double));
    // diagonal vectors, d2 is used by the solver function
    complex double* d = (complex double*) calloc(nx, sizeof(complex double));
    complex double* d2 = (complex double*) calloc(nx, sizeof(complex double));

    double ki = 0.02;
    double kf = 2.5;
    double dk = (kf - ki) / 100;
    double k_half_width = 2 / sqrt(2 * std * std);
    for (double curr_k = ki; curr_k <= kf; curr_k += dk) {
        init_x(xpos);
        init_psi0(psi, xpos, curr_k);
        Tprint(trans_data, psi, 0);

        init_V(V, xpos);
        init_d(d, d2, V);

        for (int n = 1; n < 324; n++) {
            next_phi(psi, phi, V);
            solver(psi, phi, d, d2, 1);
            Tprint(trans_data, psi, dt * n);
        }
        Fprint(F_data, psi, curr_k, k_half_width);
        for (int n = 325; n < nt; n++) {
            next_phi(psi, phi, V);
            solver(psi, phi, d, d2, 1);
            Tprint(trans_data, psi, dt * n);
        }
        fprintf(trans_data, "end\n");
    }

    free(xpos);
    free(psi);
    free(phi);
    free(V);
    free(normed);
    free(d);
    free(d2);
    fclose(trans_data);
    fclose(F_data);
}
