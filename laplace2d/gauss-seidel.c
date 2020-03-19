#include <stdio.h>
#include <math.h>

#define N 64

double phi[N+2][N+2];
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
        // gs iteration
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                phi[i][j] = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]) / 4;
        // boundary conditions
        for (int i = 1; i <= N; i++) {
            phi[i][0] = -phi[i][1];
            phi[i][N+1] = 2 - phi[i][N];
        }
        for (int j = 1; j <= N; j++) {
            phi[0][j] = -phi[1][j];
            phi[N+1][j] = -phi[N][j];
        }

        // calculate residual
        res = 0;
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                res += fabs(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4*phi[i][j]) / (h*h);
        res /= N * N;

        if (iter == 0)
            res1 = res;
        resnorm = res / res1;
        fprintf(fp, "%e\n", resnorm);

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