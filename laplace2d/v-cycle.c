#include <stdio.h>
#include <math.h>

#define N 64

double phi[N+2][N+2], phi_2h[N/2+2][N/2+2], phi_4h[N/4+2][N/4+2];
const double h = 1. / N;

int main(void) {
    int iter = 0;
    double res, res1 = 1, resnorm = 1;

    FILE *fp = fopen("residual.csv", "w");

    // initialize
    for (int i = 1; i <= N; i++) {
        phi[i][0] = -phi[i][1];
        phi[i][N+1] = 2 - phi[i][N];
    }
    for (int j = 1; j <= N; j++) {
        phi[0][j] = -phi[1][j];
        phi[N+1][j] = -phi[N][j];
    }

    while (resnorm > 1e-6) {
        // iteration on h grid
        iter_h_grid(5);
    }

    return 0;
}

void iter_h_grid(int iter_no) {
    for (int k = 0; k < iter_no; k++) {
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                phi[i][j] = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]) / 4;
        for (int i = 1; i <= N; i++) {
            phi[i][0] = -phi[i][1];
            phi[i][N+1] = 2 - phi[i][N];
        }
        for (int j = 1; j <= N; j++) {
            phi[0][j] = -phi[1][j];
            phi[N+1][j] = -phi[N][j];
        }
    }
}