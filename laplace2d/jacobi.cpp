#include <iostream>
#include <cmath>
#include <algorithm>

constexpr int N = 1024;
constexpr double phi_l = 0, phi_r = 0, phi_t = 1, phi_b = 0;

double a[N][N], b[N][N];
double (*phi)[N] = a, (*phi_new)[N] = b;

int main(void) {
    int iter = 1;
    double num, den, res;
    double res1 = 1, resnorm = 1;

    while (resnorm > 1e-6) {
        res = 0;

        for (int i = 0; i <= N-1; i++)
            for (int j = 0; j <= N-1; j++) {
                num = 0;
                den = 4;

                num += i == N-1 ? (den++, 2*phi_r) : phi[i+1][j];
                num += i == 0 ? (den++, 2*phi_l) : phi[i-1][j];
                num += j == N-1 ? (den++, 2*phi_t) : phi[i][j+1];
                num += j == 0 ? (den++, 2*phi_b) : phi[i][j-1];

                res = std::abs(phi[i][j] - num / den);
                phi_new[i][j] = num / den;
            }

        res /= N * N;
        if (iter == 1)
            res1 = res;
        resnorm = res / res1;
        std::cout << resnorm << std::endl;

        std::swap(phi, phi_new);

        iter++;
    }

    std::cout << iter;
}