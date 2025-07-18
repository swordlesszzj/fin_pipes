#include "SymmetricTridiagonalQR.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

using namespace std;

SymmetricTridiagonalQR::SymmetricTridiagonalQR(int size) : n(size) {
    d.resize(n);
    e.resize(n - 1);
    Q.resize(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        Q[i][i] = 1.0;
    }
}

void SymmetricTridiagonalQR::setDiagonal(int i, double value) {
    if (i >= 0 && i < n) d[i] = value;
}

void SymmetricTridiagonalQR::setOffDiagonal(int i, double value) {
    if (i >= 0 && i < n - 1) e[i] = value;
}

void SymmetricTridiagonalQR::setMatrix(const vector<double>& diag, const vector<double>& offDiag) {
    if (diag.size() != n || offDiag.size() != n - 1)
        throw invalid_argument("Matrix size mismatch.");
    d = diag;
    e = offDiag;
}

void SymmetricTridiagonalQR::printMatrix() const {
    cout << "Tridiagonal Matrix A:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j)
                cout << setw(10) << d[i];
            else if (abs(i - j) == 1)
                cout << setw(10) << e[min(i, j)];
            else
                cout << setw(10) << 0.0;
        }
        cout << endl;
    }
}

void SymmetricTridiagonalQR::printMatrixInfo() const {
    cout << "Diagonal (d): ";
    for (auto val : d) cout << val << " ";
    cout << "\nOff-diagonal (e): ";
    for (auto val : e) cout << val << " ";
    cout << endl;
}

bool SymmetricTridiagonalQR::isConverged(const vector<double>& diag, const vector<double>& off) const {
    for (double val : off) {
        if (fabs(val) > EPS) return false;
    }
    return true;
}

void SymmetricTridiagonalQR::deflateMatrix(vector<double>& diag, vector<double>& off) {
    while (n > 1 && fabs(off[n - 2]) < EPS) {
        n--;
        diag.pop_back();
        off.pop_back();
    }
}

double SymmetricTridiagonalQR::computeWilkinsonShift(const vector<double>& diag, const vector<double>& off) {
    int m = diag.size() - 1;
    double d1 = diag[m - 1];
    double d2 = diag[m];
    double e1 = off[m - 1];
    double delta = (d1 - d2) / 2.0;
    double sign = delta >= 0 ? 1 : -1;
    return d2 - (e1 * e1) / (fabs(delta) + sqrt(delta * delta + e1 * e1)) * sign;
}

void SymmetricTridiagonalQR::performQRStep(vector<double>& diag, vector<double>& off) {
    int m = diag.size();
    double mu = computeWilkinsonShift(diag, off);
    double x = diag[0] - mu;
    double z = off[0];

    for (int k = 0; k < m - 1; ++k) {
        double r = hypot(x, z);
        double c = x / r;
        double s = z / r;

        double dk = diag[k];
        double dk1 = diag[k + 1];
        double ek = off[k];

        diag[k] = c * c * dk + s * s * dk1 - 2 * s * c * ek;
        diag[k + 1] = s * s * dk + c * c * dk1 + 2 * s * c * ek;
        off[k] = s * c * (dk - dk1) + (c * c - s * s) * ek;

        if (k < m - 2) {
            double ek1 = off[k + 1];
            x = off[k];
            z = -s * ek1;
            off[k + 1] = c * ek1;
        }
    }
}

void SymmetricTridiagonalQR::performQRStepWithVectors(vector<double>& diag, vector<double>& off) {
    int m = diag.size();
    double mu = computeWilkinsonShift(diag, off);
    double x = diag[0] - mu;
    double z = off[0];

    for (int k = 0; k < m - 1; ++k) {
        double r = hypot(x, z);
        double c = x / r;
        double s = z / r;

        double dk = diag[k];
        double dk1 = diag[k + 1];
        double ek = off[k];

        diag[k] = c * c * dk + s * s * dk1 - 2 * s * c * ek;
        diag[k + 1] = s * s * dk + c * c * dk1 + 2 * s * c * ek;
        off[k] = s * c * (dk - dk1) + (c * c - s * s) * ek;

        if (k < m - 2) {
            double ek1 = off[k + 1];
            x = off[k];
            z = -s * ek1;
            off[k + 1] = c * ek1;
        }

        // Apply rotation to eigenvector matrix
        for (int i = 0; i < n; ++i) {
            double temp = c * Q[i][k] + s * Q[i][k + 1];
            Q[i][k + 1] = -s * Q[i][k] + c * Q[i][k + 1];
            Q[i][k] = temp;
        }
    }
}

std::vector<double> SymmetricTridiagonalQR::computeEigenvalues() {
    std::vector<double> diag = d;
    std::vector<double> off = e;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        if (isConverged(diag, off)) break;
        deflateMatrix(diag, off);
        performQRStep(diag, off);
    }

    return diag;
}

std::pair<std::vector<double>, std::vector<std::vector<double>>> SymmetricTridiagonalQR::computeEigenvaluesAndVectors() {
    std::vector<double> diag = d;
    std::vector<double> off = e;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        if (isConverged(diag, off)) break;
        deflateMatrix(diag, off);
        performQRStepWithVectors(diag, off);
    }

    return {diag, Q};
}

double SymmetricTridiagonalQR::estimateConditionNumber(const std::vector<double>& eigenvalues) const {
    double maxVal = *std::max_element(eigenvalues.begin(), eigenvalues.end());
    double minVal = *std::min_element(eigenvalues.begin(), eigenvalues.end());
    return fabs(maxVal / minVal);
}

double SymmetricTridiagonalQR::computeCharacteristicPolynomial(double lambda) const {
    if (n == 1) return d[0] - lambda;
    std::vector<double> p(n + 1);
    p[0] = 1.0;
    p[1] = d[0] - lambda;
    for (int i = 2; i <= n; i++) {
        p[i] = (d[i - 1] - lambda) * p[i - 1] - e[i - 2] * e[i - 2] * p[i - 2];
    }
    return p[n];
}

void SymmetricTridiagonalQR::verifyEigenvalues(const std::vector<double>& eigenvalues) const {
    cout << "Verification (Characteristic Polynomial):\n";
    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        double lambda = eigenvalues[i];
        double val = computeCharacteristicPolynomial(lambda);
        cout << "λ_" << i + 1 << " = " << lambda << ", P(λ) = " << val << endl;
    }
}
