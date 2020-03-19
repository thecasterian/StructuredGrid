#include <stdio.h>
#include <math.h>

#define N 64

double a[N+2][N+2], b[N+2][N+2];
double (*phi)[N+2] = a, (*phi_new)[N+2] = b;
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
        // jacobi iteration
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                phi_new[i][j] = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]) / 4;
        // boundary conditions
        for (int i = 1; i <= N; i++) {
            phi_new[i][0] = -phi_new[i][1];
            phi_new[i][N+1] = 2 - phi_new[i][N];
        }
        for (int j = 1; j <= N; j++) {
            phi_new[0][j] = -phi_new[1][j];
            phi_new[N+1][j] = -phi_new[N][j];
        }

        // calculate residual
        res = 0;
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                res += fabs(phi_new[i+1][j]+phi_new[i-1][j]+phi_new[i][j+1]+phi_new[i][j-1]-4*phi_new[i][j]) / (h*h);
        res /= N * N;

        if (iter == 1)
            res1 = res;
        resnorm = res / res1;
        fprintf(fp, "%e\n", resnorm);

        double (*tmp)[N+2] = phi;
        phi = phi_new;
        phi_new = tmp;

        iter++;
    }

    printf("%d\n", iter);
    fclose(fp);

    fp = fopen("phi.csv", "w");
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            fprintf(fp, "%e", phi[i][j]);
            if (j < N)
                fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}