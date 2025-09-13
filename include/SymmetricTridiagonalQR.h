#ifndef SYMMETRICTRIDIAGONALQR_H
#define SYMMETRICTRIDIAGONALQR_H

#include <vector>
#include <iostream>
#include <iomanip>  
#include <algorithm>
#include <cmath>
#include <stdexcept>

class SymmetricTridiagonalQR {
private:
    std::vector<double> d;  // 主对角线元素
    std::vector<double> e;  // 次对角线元素（上下对角线相同）
    std::vector<std::vector<double>> Q;  // 特征向量矩阵
    int n;

    const double EPS = 1e-200;
    const int MAX_ITER = 11000;

public:
    SymmetricTridiagonalQR(int size);

    void setDiagonal(int i, double value);
    void setOffDiagonal(int i, double value);
    void setMatrix(const std::vector<double>& diagonal, const std::vector<double>& offDiagonal);
    void printMatrix() const;

    std::pair<std::vector<double>, std::vector<std::vector<double>>> computeEigenvaluesAndVectors();
    std::vector<double> computeEigenvalues();

    void printMatrixInfo() const;
    void verifyEigenvalues(const std::vector<double>& eigenvalues) const;

private:
    double computeWilkinsonShift(const std::vector<double>& diag, const std::vector<double>& off);
    void performQRStep(std::vector<double>& diag, std::vector<double>& off);
    void performQRStepWithVectors(std::vector<double>& diag, std::vector<double>& off);
    bool isConverged(const std::vector<double>& diag, const std::vector<double>& off) const;
    void deflateMatrix(std::vector<double>& diag, std::vector<double>& off);
    double estimateConditionNumber(const std::vector<double>& eigenvalues) const;
    double computeCharacteristicPolynomial(double lambda) const;
};

#endif // SYMMETRICTRIDIAGONALQR_H
