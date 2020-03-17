#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 64

double phi[N+2][N+2];
double res_h[N+2][N+2], res_2h[N/2+2][N/2+2], res_4h[N/4+2][N/4+2];
double err_h[N+2][N+2], err_2h[N/2+2][N/2+2], err_4h[N/4+2][N/4+2];
const double h = 1. / N;

void smoothing(int iter_no);
void calc_res_h(void);
void rest_to_2h(void);
void calc_err_2h(int iter_no);
void prol_to_h(void);
void correct_h(void);

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
        smoothing(5);
        // calculate residual on h grid
        calc_res_h();
        // restrict residual to 2h grid
        rest_to_2h();
        // calculate error vector on 2h grid
        calc_err_2h(10);
        // prolongate error vector to h grid
        prol_to_h();
        // correct on h grid
        correct_h();
        // iteration on h grid
        smoothing(5);

        // calculate residual
        res = 0;
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                res += fabs(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4*phi[i][j]) / (h*h);
        res /= N * N;

        if (iter == 0)
            res1 = res;
        resnorm = res / res1;
        printf("%e\n", resnorm);

        iter++;
    }

    printf("%d\n", iter);
    fclose(fp);
    return 0;
}

void smoothing(int iter_no) {
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

void calc_res_h(void) {
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            res_h[i][j] = -(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4*phi[i][j]);
}

void rest_to_2h(void) {
    for (int i = 1; i <= N/2; i++)
        for (int j = 1; j <= N/2; j++)
            res_2h[i][j] = (res_h[2*i-1][2*j-1]+res_h[2*i][2*j-1]+res_h[2*i-1][2*j]+res_h[2*i][2*j]) / 4;
}

void calc_err_2h(int iter_no) {
    // initialize err_2h to zero
    memset(err_2h, 0, sizeof(double)*(N/2+2)*(N/2+2));
    // gs iteration
    for (int k = 0; k < iter_no; k++) {
        for (int i = 1; i <= N/2; i++)
            for (int j = 1; j <= N/2; j++)
                err_2h[i][j] = (err_2h[i+1][j]+err_2h[i-1][j]+err_2h[i][j+1]+err_2h[i][j-1]) / 4;
        for (int i = 1; i <= N/2; i++) {
            err_2h[i][0] = -err_2h[i][1];
            err_2h[i][N/2+1] = -err_2h[i][N/2];
        }
        for (int j = 0; j <= N/2+1; j++) {
            err_2h[0][j] = -err_2h[1][j];
            err_2h[N/2+1][j] = -err_2h[N/2][j];
        }
    }
}

void prol_to_h(void) {
    for (int i = 1; 2*i <= N; i++)
        for (int j = 1; 2*j <= N; j++)
            err_h[2*i][2*j] = (9*err_2h[i][j] + 3*err_2h[i+1][j] + 3*err_2h[i][j+1] + err_2h[i+1][j+1]) / 16;
    for (int i = 0; 2*i+1 <= N; i++)
        for (int j = 1; 2*j <= N; j++)
            err_h[2*i+1][2*j] = (3*err_2h[i][j] + 9*err_2h[i+1][j] + 3*err_2h[i][j+1] + err_2h[i+1][j+1]) / 16;
    for (int i = 1; 2*i <= N; i++)
        for (int j = 0; 2*j+1 <= N; j++)
            err_h[2*i][2*j+1] = (3*err_2h[i][j] + err_2h[i+1][j] + 9*err_2h[i][j+1] + 3*err_2h[i+1][j+1]) / 16;
    for (int i = 0; 2*i+1 <= N; i++)
        for (int j = 0; 2*j+1 <= N; j++)
            err_h[2*i+1][2*j+1] = (err_2h[i][j] + 3*err_2h[i+1][j] + 3*err_2h[i][j+1] + 9*err_2h[i+1][j+1]) / 16;
}

void correct_h(void) {
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            phi[i][j] += err_h[i][j];
    for (int i = 1; i <= N; i++) {
        phi[i][0] = -phi[i][1];
        phi[i][N+1] = 2 - phi[i][N];
    }
    for (int j = 1; j <= N; j++) {
        phi[0][j] = -phi[1][j];
        phi[N+1][j] = -phi[N][j];
    }
}