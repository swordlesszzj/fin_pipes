#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <utility>
#include "SymmetricTridiagonalQR.h" 

using namespace std;

class SymmetricTridiagonalQR; // 前置声明，外部类需自行包含头文件

class Matrix {
public:
    int rows;
    int cols;
    vector<vector<double>> data;

    bool isSquare() const;

    // 构造函数
    Matrix(int r, int c);
    Matrix(const vector<vector<double>>& vec);

    // 运算符重载
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    friend Matrix operator*(double scalar, const Matrix& m);
    Matrix operator/(double scalar) const;

    // 矩阵访问操作符
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    // 打印函数
    friend ostream& operator<<(ostream& os, const Matrix& m);

    Matrix transpose() const;

    // 静态函数 解线性方程 Ax = B
    static Matrix Solve_Matrix(const Matrix& A, const Matrix& B);

    // Lanczos迭代法
    pair<vector<double>, Matrix> lanczosIteration(const Matrix& K,const Matrix& M,int q,double tol = 1e-12);

    // Cholesky分解
    Matrix choleskyDecomposition(const Matrix& M);

    // 子空间迭代法
    pair<vector<double>, Matrix> subspaceIteration(const Matrix& K, const Matrix& M, int q, double tol = 1e-6, int max_iter = 100);

    // 列向量获取和设置
    Matrix getColumn(int col) const;
    void setColumn(int col, const Matrix& vec);

    // 矩阵信息获取
    int getRows() const;
    int getCols() const;

    // 其它成员函数（如 StandardEigenSolver）如需添加请声明
};
