#include <stdlib.h>
#include "prop.h"

// prints x position, norm of psi, and psi parts
static void write_vec(FILE* data, const double* normed,
        const complex double* psi, const double* xpos) {
    for (int j = 0; j < nx; j++) {
        fprintf(data, "%lf %.6lf %.6lf %.6lf\n", xpos[j], normed[j],
                creal(psi[j]), cimag(psi[j]));
    }
    fprintf(data, "end\n");
}

int main() {
    FILE* data = fopen(".used/datafiles/wave.dat", "w");
    FILE* potential = fopen(".used/datafiles/potential.dat", "w");
    double* xpos = (double*) calloc(nx, sizeof(double));
    complex double* psi = (complex double*) calloc(nx, sizeof(complex double));
    double* normed = (double*) calloc(nx, sizeof(double));
    complex double* phi = (complex double*) calloc(nx, sizeof(complex double));
    double* V = (double*) calloc(nx, sizeof(double));
    // diagonal vectors, d2 is used by the solver function
    complex double* d = (complex double*) calloc(nx, sizeof(complex double));
    complex double* d2 = (complex double*) calloc(nx, sizeof(complex double));

    float k, dt;
    scanf("k=%f,dt=%f", &k, &dt);
    init_x(xpos);
    init_psi0(psi, xpos, k);
    norm(normed, psi);
    write_vec(data, normed, psi, xpos);

    init_V(V, xpos);
    init_d(d, d2, V, dt);
    for (int j = 0; j < nx; j++) {
        fprintf(potential, "%lf %lf\n", xpos[j], V[j]);
    }

    for (int n = 1; n < nt; n++) {
        next_phi(psi, phi, V, dt);
        solver(psi, phi, d, d2, 1, dt);
        norm(normed, psi);
        write_vec(data, normed, psi, xpos);
    }

    free(xpos);
    free(psi);
    free(phi);
    free(V);
    free(normed);
    free(d);
    free(d2);
    fclose(data);
    fclose(potential);
}
