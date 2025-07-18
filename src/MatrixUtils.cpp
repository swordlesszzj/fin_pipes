#include "MatrixUtils.h"

// 对称矩阵补充上半，a为12阶矩阵，b=12
void symmetry(double a[][12], double b) {
    for (int i = 0; i <= b - 2; i++) {
        for (int j = i + 1; j <= b - 1; j++) {
            a[i][j] = a[j][i];
        }
    }
}

// 矩阵乘法，结果存MultiResult，均为12阶矩阵
void Multi(double a[][12], double b[][12], double c[][12], double MultiResult[][12]) {
    double MultiResult1[12][12]{};
    std::memset(MultiResult, 0, sizeof(double) * 12 * 12);
    std::memset(MultiResult1, 0, sizeof(double) * 12 * 12);

    for (int i = 0; i <= 11; i++) {
        for (int j = 0; j <= 11; j++) {
            for (int k = 0; k <= 11; k++) {
                MultiResult1[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    for (int i = 0; i <= 11; i++) {
        for (int j = 0; j <= 11; j++) {
            for (int k = 0; k <= 11; k++) {
                MultiResult[i][j] += MultiResult1[i][k] * c[k][j];
            }
        }
    }
}
