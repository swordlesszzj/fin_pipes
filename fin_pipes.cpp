#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <iomanip>
#include <string>
#include <fstream>
#include "./json.hpp"
using namespace std;
using json = nlohmann::json;
const float Pi = 3.1415926;

double MultiResult[12][12];



struct Pipe {

    int Id;
    int StartNode;
    int EndNode;
    float DeltaX;
    float DeltaY;
    float DeltaZ;
    float Diameter;//外径
    float WallThickness;
    float ElasticModulus;
    float PoissionRatio;
    float Density;

    //生成中间变量，读入数据时不需要管
    float I_x;//转动惯量
    float J_x, J_y;//惯性矩
    float ShearModulus;
    float Mass;
    float Length;
    float Aera;
};
double K1[12][12], M1[12][12], K[12][12], M[12][12], c[3][3], C[12][12],CT[12][12];//旋转前后矩阵，旋转因子、旋转矩阵、转置后的旋转矩阵

void symmetry(double a[][12],double b) {//对称矩阵补充上半,参数是指针和矩阵阶数
    for (int i = 0; i <= b - 2; i++ ){
        for (int j = i + 1; j <= b - 1; j++) {
            a[i][j] = a[j][i];
        }
    }
}


void Multi(double a[][12], double b[][12], double c[][12]) {
    memset(MultiResult, 0, sizeof(MultiResult));
    double MultiResult1[12][12];
    memset(MultiResult1, 0, sizeof(MultiResult1));
    for (int i = 0; i <= 11; i++) {
        for (int j = 0; j <= 11; j++) {
            for (int k = 0; k <= 11; k++) {
                MultiResult1[i][j] += a[i][k] * b[k][j];
            }
        }
    }
   // symmetry(MultiResult1, 12);
    for (int i = 0; i <= 11; i++) {
        for (int j = 0; j <= 11; j++) {
            for (int k = 0; k <= 11; k++) {
                MultiResult[i][j] += MultiResult1[i][k] * c[k][j];
            }
        }
    }
    //symmetry(MultiResult, 12);
}

class SymmetricTridiagonalQR {
private:
    std::vector<double> d;  // 主对角线元素
    std::vector<double> e;  // 次对角线元素（上下对角线相同）
    std::vector<std::vector<double>> Q;  // 特征向量矩阵
    int n;

    const double EPS = 1e-15;
    const int MAX_ITER = 1000;

public:
    SymmetricTridiagonalQR(int size) : n(size) {
        d.resize(n);
        e.resize(n - 1);
        Q.resize(n, std::vector<double>(n, 0.0));
        // 初始化Q为单位矩阵
        for (int i = 0; i < n; i++) {
            Q[i][i] = 1.0;
        }
    }

    // 设置主对角线元素
    void setDiagonal(int i, double value) {
        if (i >= 0 && i < n) d[i] = value;
    }

    // 设置次对角线元素（对称矩阵，上下对角线相同）
    void setOffDiagonal(int i, double value) {
        if (i >= 0 && i < n - 1) e[i] = value;
    }

    // 批量设置矩阵
    void setMatrix(const std::vector<double>& diagonal, const std::vector<double>& offDiagonal) {
        if (diagonal.size() != n || offDiagonal.size() != n - 1) {
            throw std::invalid_argument("数组大小不匹配");
        }
        d = diagonal;
        e = offDiagonal;
    }

    // 打印矩阵
    void printMatrix() const {
        std::cout << "对称三对角矩阵：\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    std::cout << std::setw(10) << std::fixed << std::setprecision(4) << d[i];
                }
                else if (i == j + 1) {
                    std::cout << std::setw(10) << std::fixed << std::setprecision(4) << e[j];
                }
                else if (i == j - 1) {
                    std::cout << std::setw(10) << std::fixed << std::setprecision(4) << e[i];
                }
                else {
                    std::cout << std::setw(10) << std::fixed << std::setprecision(4) << 0.0;
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // 主QR算法求解特征值和特征向量
    std::pair<std::vector<double>, std::vector<std::vector<double>>> computeEigenvaluesAndVectors() {
        std::vector<double> diag = d;
        std::vector<double> off = e;

        // 重新初始化Q为单位矩阵
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        std::cout << "开始QR迭代...\n";

        for (int iter = 0; iter < MAX_ITER; iter++) {
            // 检查收敛性
            if (isConverged(diag, off)) {
                std::cout << "QR算法在第 " << iter << " 次迭代后收敛\n";
                break;
            }

            // 计算Wilkinson位移
            double shift = computeWilkinsonShift(diag, off);

            // 应用位移
            for (int i = 0; i < n; i++) {
                diag[i] -= shift;
            }

            // 执行QR步骤并更新特征向量
            performQRStepWithVectors(diag, off);

            // 恢复位移
            for (int i = 0; i < n; i++) {
                diag[i] += shift;
            }

            // 检查并处理收敛的特征值
            deflateMatrix(diag, off);

            // 每100次迭代输出一次进度
            if ((iter + 1) % 100 == 0) {
                std::cout << "完成第 " << iter + 1 << " 次迭代\n";
            }
        }

        // 对特征值和特征向量进行排序
        std::vector<int> indices(n);
        for (int i = 0; i < n; i++) indices[i] = i;

        std::sort(indices.begin(), indices.end(), [&](int i, int j) {
            return diag[i] < diag[j];
            });

        std::vector<double> sorted_eigenvalues(n);
        std::vector<std::vector<double>> sorted_eigenvectors(n, std::vector<double>(n));

        for (int i = 0; i < n; i++) {
            sorted_eigenvalues[i] = diag[indices[i]];
            for (int j = 0; j < n; j++) {
                sorted_eigenvectors[j][i] = Q[j][indices[i]];
            }
        }

        return { sorted_eigenvalues, sorted_eigenvectors };
    }

    // 仅计算特征值（保持原有功能）
    std::vector<double> computeEigenvalues() {
        std::vector<double> diag = d;
        std::vector<double> off = e;

        std::cout << "开始QR迭代...\n";

        for (int iter = 0; iter < MAX_ITER; iter++) {
            // 检查收敛性
            if (isConverged(diag, off)) {
                std::cout << "QR算法在第 " << iter << " 次迭代后收敛\n";
                break;
            }

            // 计算Wilkinson位移
            double shift = computeWilkinsonShift(diag, off);

            // 应用位移
            for (int i = 0; i < n; i++) {
                diag[i] -= shift;
            }

            // 执行QR步骤
            performQRStep(diag, off);

            // 恢复位移
            for (int i = 0; i < n; i++) {
                diag[i] += shift;
            }

            // 检查并处理收敛的特征值
            deflateMatrix(diag, off);

            // 每100次迭代输出一次进度
            if ((iter + 1) % 100 == 0) {
                std::cout << "完成第 " << iter + 1 << " 次迭代\n";
            }
        }

        // 排序特征值
        std::sort(diag.begin(), diag.end());
        return diag;
    }

private:
    // 计算Wilkinson位移（对称三对角矩阵）
    double computeWilkinsonShift(const std::vector<double>& diag, const std::vector<double>& off) {
        if (n == 1) return diag[0];

        // 使用右下角2x2子矩阵
        double a = diag[n - 2];
        double b = diag[n - 1];
        double c = off[n - 2];

        // 计算2x2对称矩阵 [a c; c b] 的特征值
        double trace = a + b;
        double det = a * b - c * c;
        double discriminant = trace * trace - 4 * det;

        if (discriminant < 0) {
            return b;
        }

        double sqrt_disc = std::sqrt(discriminant);
        double lambda1 = (trace + sqrt_disc) / 2.0;
        double lambda2 = (trace - sqrt_disc) / 2.0;

        // 选择离b更近的特征值作为位移
        return (std::abs(lambda1 - b) < std::abs(lambda2 - b)) ? lambda1 : lambda2;
    }

    // 执行一步QR分解和更新（仅特征值）
    void performQRStep(std::vector<double>& diag, std::vector<double>& off) {
        // 使用Givens旋转进行QR分解
        std::vector<double> c(n - 1), s(n - 1);

        // 前向消元：将矩阵变成上三角
        for (int i = 0; i < n - 1; i++) {
            double a = diag[i];
            double b = off[i];

            // 计算Givens旋转参数
            double r = std::sqrt(a * a + b * b);

            if (std::abs(r) < EPS) {
                c[i] = 1.0;
                s[i] = 0.0;
                continue;
            }

            c[i] = a / r;
            s[i] = b / r;

            // 更新矩阵元素
            diag[i] = r;

            // 更新下一行
            if (i < n - 1) {
                double temp1 = c[i] * diag[i + 1] + s[i] * (i < n - 2 ? off[i + 1] : 0.0);
                double temp2 = -s[i] * diag[i + 1] + c[i] * (i < n - 2 ? off[i + 1] : 0.0);

                diag[i + 1] = temp1;
                if (i < n - 2) {
                    off[i + 1] = temp2;
                }
            }
        }

        // 后向更新：RQ相乘恢复三对角形式
        for (int i = n - 2; i >= 0; i--) {
            double temp_d = c[i] * diag[i] + s[i] * off[i];
            double temp_e = -s[i] * diag[i] + c[i] * off[i];

            diag[i] = temp_d;
            off[i] = temp_e;

            if (i > 0) {
                double temp = c[i] * diag[i - 1] + s[i] * off[i - 1];
                off[i - 1] = -s[i] * diag[i - 1] + c[i] * off[i - 1];
                diag[i - 1] = temp;
            }
        }
    }

    // 执行一步QR分解和更新（包含特征向量）
    void performQRStepWithVectors(std::vector<double>& diag, std::vector<double>& off) {
        std::vector<double> c(n - 1), s(n - 1);

        // 前向消元：将矩阵变成上三角
        for (int i = 0; i < n - 1; i++) {
            double a = diag[i];
            double b = off[i];

            double r = std::sqrt(a * a + b * b);

            if (std::abs(r) < EPS) {
                c[i] = 1.0;
                s[i] = 0.0;
                continue;
            }

            c[i] = a / r;
            s[i] = b / r;

            diag[i] = r;

            if (i < n - 1) {
                double temp1 = c[i] * diag[i + 1] + s[i] * (i < n - 2 ? off[i + 1] : 0.0);
                double temp2 = -s[i] * diag[i + 1] + c[i] * (i < n - 2 ? off[i + 1] : 0.0);

                diag[i + 1] = temp1;
                if (i < n - 2) {
                    off[i + 1] = temp2;
                }
            }
        }

        // 后向更新：RQ相乘恢复三对角形式
        for (int i = n - 2; i >= 0; i--) {
            double temp_d = c[i] * diag[i] + s[i] * off[i];
            double temp_e = -s[i] * diag[i] + c[i] * off[i];

            diag[i] = temp_d;
            off[i] = temp_e;

            if (i > 0) {
                double temp = c[i] * diag[i - 1] + s[i] * off[i - 1];
                off[i - 1] = -s[i] * diag[i - 1] + c[i] * off[i - 1];
                diag[i - 1] = temp;
            }

            // 更新特征向量矩阵
            for (int j = 0; j < n; j++) {
                double temp1 = c[i] * Q[j][i] + s[i] * Q[j][i + 1];
                double temp2 = -s[i] * Q[j][i] + c[i] * Q[j][i + 1];
                Q[j][i] = temp1;
                Q[j][i + 1] = temp2;
            }
        }
    }

    // 检查收敛性
    bool isConverged(const std::vector<double>& diag, const std::vector<double>& off) const {
        for (int i = 0; i < n - 1; i++) {
            double tolerance = EPS * (std::abs(diag[i]) + std::abs(diag[i + 1]));
            if (std::abs(off[i]) > tolerance) {
                return false;
            }
        }
        return true;
    }

    // 处理收敛的特征值（矩阵分解）
    void deflateMatrix(std::vector<double>& diag, std::vector<double>& off) {
        for (int i = 0; i < n - 1; i++) {
            double tolerance = EPS * (std::abs(diag[i]) + std::abs(diag[i + 1]));
            if (std::abs(off[i]) <= tolerance) {
                off[i] = 0.0;
            }
        }
    }

    // 计算条件数估计
    double estimateConditionNumber(const std::vector<double>& eigenvalues) const {
        if (eigenvalues.empty()) return 1.0;

        double max_val = *std::max_element(eigenvalues.begin(), eigenvalues.end());
        double min_val = *std::min_element(eigenvalues.begin(), eigenvalues.end());

        if (std::abs(min_val) < EPS) return 1e16;

        return std::abs(max_val / min_val);
    }

public:
    // 获取矩阵的基本信息
    void printMatrixInfo() const {
        std::cout << "矩阵维度: " << n << "x" << n << "\n";
        std::cout << "主对角线元素: ";
        for (int i = 0; i < n; i++) {
            std::cout << std::fixed << std::setprecision(4) << d[i];
            if (i < n - 1) std::cout << ", ";
        }
        std::cout << "\n次对角线元素: ";
        for (int i = 0; i < n - 1; i++) {
            std::cout << std::fixed << std::setprecision(4) << e[i];
            if (i < n - 2) std::cout << ", ";
        }
        std::cout << "\n\n";
    }

    // 验证特征值精度
    void verifyEigenvalues(const std::vector<double>& eigenvalues) const {

        for (int i = 0; i < eigenvalues.size(); i++) {
            double lambda = eigenvalues[i];
            double det = computeCharacteristicPolynomial(lambda);

        }

        double cond = estimateConditionNumber(eigenvalues);

    }

private:
    // 计算特征多项式的值 det(A - λI)
    double computeCharacteristicPolynomial(double lambda) const {
        if (n == 1) return d[0] - lambda;

        std::vector<double> p(n + 1);
        p[0] = 1.0;
        p[1] = d[0] - lambda;

        for (int i = 2; i <= n; i++) {
            p[i] = (d[i - 1] - lambda) * p[i - 1] - e[i - 2] * e[i - 2] * p[i - 2];
        }

        return p[n];
    }
};

class Matrix {
public:
    int rows;
    int cols;
    vector<vector<double>> data;
    bool isSquare() const {
        return rows == cols;
    }
    // 构造函数
    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(rows, vector<double>(cols, 0.0));
    }

    Matrix(const vector<vector<double>>& vec) {
        if (vec.empty() || vec[0].empty()) {
            throw invalid_argument("无效矩阵数据");
        }
        rows = vec.size();
        cols = vec[0].size();
        for (const auto& row : vec) {
            if (row.size() != cols) {
                throw invalid_argument("矩阵行长度不一致");
            }
        }
        data = vec;
    }

    // 运算符重载
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("矩阵维度不匹配");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("矩阵维度不匹配");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw invalid_argument("矩阵维度不匹配");
        }
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < cols; ++k) {
                    sum += data[i][k] * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }

    friend Matrix operator*(double scalar, const Matrix& m) {
        return m * scalar;
    }

    Matrix operator/(double scalar) const {
        if (scalar == 0) {
            throw invalid_argument("除数不能为零");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] / scalar;
            }
        }
        return result;
    }

    
    // 矩阵访问操作符
    double& operator()(int i, int j) {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw out_of_range("索引越界");
        }
        return data[i][j];
    }

    const double& operator()(int i, int j) const {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw out_of_range("索引越界");
        }
        return data[i][j];
    }

    // 打印函数
    friend ostream& operator<<(ostream& os, const Matrix& m) {
        for (int i = 0; i < m.rows; ++i) {
            for (int j = 0; j < m.cols; ++j) {
                os << m.data[i][j];
                if (j < m.cols - 1) os << "\t";
            }
            os << endl;
        }
        return os;
    }

    Matrix transpose() const {
        vector<vector<double>> transposed(cols, vector<double>(rows));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                transposed[j][i] = data[i][j];
            }
        }
        Matrix temp(transposed);
        return temp;
    }
// 解线性方程 Ax = B
    static Matrix Solve_Matrix(const Matrix& A, const Matrix& B) {
        // 检查输入有效性
        if (!A.isSquare()) {
            throw invalid_argument("系数矩阵A必须为方阵");
        }
        if (A.getRows() != B.getRows()) {
            throw invalid_argument("矩阵A和B行数必须匹配");
        }

        int n = A.getRows();
        int m = B.getCols();
        Matrix x(n, m);

        // 复制矩阵数据作增广矩阵 
        vector<vector<double>> Ab(n, vector<double>(n + m));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ab[i][j] = A(i, j);
            }
            for (int j = 0; j < m; ++j) {
                Ab[i][n + j] = B(i, j);
            }
        }

        // 高斯消元（带部分选主元）
        const double EPS = 1e-8;
        for (int k = 0; k < n; ++k) {
            // 部分选主元
            int maxRow = k;
            for (int i = k; i < n; ++i) {
                if (abs(Ab[i][k]) > abs(Ab[maxRow][k])) { maxRow = i; }
            }
            if (abs(Ab[maxRow][k]) < EPS) {
                throw domain_error("矩阵奇异，无法求解");
            }

            // 交换行
            swap(Ab[k], Ab[maxRow]);

            // 归一化主元行
            double pivot = Ab[k][k];
            for (int j = k; j < n + m; ++j) {
                Ab[k][j] /= pivot;
            }

            // 消元
            for (int i = 0; i < n; ++i) {
                if (i != k && abs(Ab[i][k]) > EPS) {
                    double factor = Ab[i][k];
                    for (int j = k; j < n + m; ++j) {
                        Ab[i][j] -= factor * Ab[k][j];
                    }
                }
            }
        }

        // 回代求解
        for (int col = 0; col < m; ++col) {
            for (int i = 0; i < n; ++i) {
                x(i, col) = Ab[i][n + col];
                for (int j = 0; j < i; ++j) {
                    x(i, col) -= Ab[i][j] * x(j, col);
                }
            }
        }

        return x;
    }
    
    pair<vector<double>, Matrix> lanczosIteration(const Matrix& K,const Matrix& M,int q,double tol = 1e-12)
    {
        int n = K.rows; 
        int r = min(2 * q, q + 8);  

        Matrix x(n, 1);  // 初始全1向量
        for (int i = 0; i < n; ++i) {
            x(i, 0) = 1.0;  
        }

        // 正则化
        Matrix temp= x.transpose() * M * x; 
        double beta = sqrt(temp(0,0));
       
        x = x / beta;

        // 存储 Lanczos 向量
        Matrix X(n, r);  // Lanczos 向量矩阵
        X.setColumn(0, x);  // 第一列

        vector<double> alpha(r, 0.0); 
        vector<double> beta_vals(r , 0.0);
        beta_vals[0] = beta;

       // Matrix x_prev(n, 1);  // 前一个向量

        //  Lanczos 迭代
        for (int i = 1; i < r; ++i) {
            // 求解 K*x_hat = M*x_{i-1}
            Matrix b = M * X.getColumn(i - 1);
            Matrix x_hat = Solve_Matrix(K, b);

            // 计算 α_{i-1} = x_hatᵀ * M * x_{i-1}
            temp= x_hat.transpose() * M * X.getColumn(i - 1);
            alpha[i - 1] = temp(0, 0);
            

            //  正交化: x_hat = x_hat - α_{i-1}x_{i-1} - β_{i-1}x_{i-2}
            if (i == 1) {
                x_hat = x_hat - alpha[i - 1] * X.getColumn(i - 1);
            }
            else {
                x_hat = x_hat - alpha[i - 1] * X.getColumn(i - 1)
                    - beta_vals[i - 1] * X.getColumn(i - 2);
            }
            /*
            // 重正交化 ???
            for (int j = 0; j < i; ++j) {
                Matrix xj = X.getColumn(j);
                Matrix Mxj = M * xj;
                double dot = 0.0;
                for (int k = 0; k < n; ++k) {
                    dot += x_hat(k, 0) * Mxj(k, 0);
                }
                x_hat = x_hat - dot * xj;
            }
            */
            // 正则化: β_i = sqrt(x_hatᵀMx_hat)
            temp = x_hat.transpose() * M * x_hat;
            beta_vals[i] = sqrt(temp(0,0));

            if (beta_vals[i] < tol) {
                // 提前终止 (剩余维度不足)
                r = i;
                break;
            }

            
            x_hat = x_hat / beta_vals[i];

            // 存储当前向量
            X.setColumn(i, x_hat); 
        }

        // 构建三对角矩阵 T
        SymmetricTridiagonalQR T(r);
        for (int i = 0; i < r; ++i) {
            T.setDiagonal(i, alpha[i]);
            if (i < r - 1) {
                T.setOffDiagonal(i,beta_vals[i + 1]);
            }
        }

        auto result = T.computeEigenvaluesAndVectors();
        auto eigenvalues = result.first;
        auto Z = result.second;

        //  求解三对角矩阵特征值问题
        //auto [eigenvalues, Z] = StandardEigenSolver(T);
        Matrix Phi = X * Z;

        vector<double> omega(r);
        for (int i = 0; i < r; ++i) {
            omega[i] = 1.0 / sqrt(eigenvalues[i]);
        }

        // 返回前q个特征对
        vector<double> omega_q(q);
        Matrix Phi_q(n, q);

        // 选取最小特征值对应的特征向量
        vector<size_t> indices(r);
        for (size_t i = 0; i < r; ++i) indices[i] = i;
        sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
            return omega[a] < omega[b];
            });

        for (int i = 0; i < q; ++i) {
            omega_q[i] = omega[indices[i]];
            Phi_q.setColumn(i, Phi.getColumn(indices[i]));
        }

        return make_pair(omega_q, Phi_q);
    }
    Matrix choleskyDecomposition(const Matrix& M) {
        int n = M.rows;
        Matrix L(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = 0.0;
                for (int k = 0; k < j; ++k) {
                    sum += L(i, k) * L(j, k);
                }
                if (i == j) {
                    L(i, j) = sqrt(M(i, i) - sum);
                }
                else {
                    L(i, j) = (M(i, j) - sum) / L(j, j);
                }
            }
        }
        return L;
    }
    pair<vector<double>, Matrix> subspaceIteration(const Matrix& K, const Matrix& M, int q, double tol = 1e-6, int max_iter = 100) {
        int n = K.rows; 
        int r = min(2 * q, q + 8); 

        Matrix X0(n, r);
        for (int i = 0; i < n; ++i) {
            X0(i, 0) = 1.0;  // 第一列全1
        }
        for (int i = 1; i < r; ++i) {
            if (i - 1 < n) X0(i - 1, i) = 1.0;  
        }

        Matrix Xk = X0;  
        vector<double> prev_eigenvalues(q, 0.0);  // 存储前一次迭代特征值
        Matrix eigenvectors(n, q);  // 最终特征向量
        
        for (int iter = 0; iter < max_iter; ++iter) {
           cout<< "迭代次数: " << iter + 1 << endl;
           cout<<M<<endl;
           cout<<Xk<<endl;
            Matrix Y = M * Xk;

            Matrix Xk1 = Solve_Matrix(K, Y);

            Matrix K_star = Xk1.transpose() * Y;

            Matrix Y1 = M * Xk1;

            Matrix M_star = Xk1.transpose() * Y1;

            Matrix L= choleskyDecomposition(M_star);
            //Matrix A = generalizedToStandard(K_star, M_star, L);

            auto result=lanczosIteration(K_star, M_star,r);
            auto eigenvalues = result.first;
            auto Z = result.second;

            //auto [eigenvalues, Z] = 

            Matrix Phi = Solve_Matrix(L.transpose(), Z);
            //取前q个：
            vector<double> current_eigenvalues(q);
            Matrix Phi_q(r, q);
            for (int i = 0; i < q; ++i) {
                current_eigenvalues[i] = eigenvalues[i];
                for (int j = 0; j < r; ++j) {
                    Phi_q(j, i) = Phi(j, i);
                }
            }

            //验敛
            bool converged = true;
            for (int i = 0; i < q; ++i) {
                if (abs(current_eigenvalues[i] - prev_eigenvalues[i]) > tol * abs(current_eigenvalues[i])) {
                    converged = false;
                    break;
                }
            }

            if (converged) {
                // 计算最终特征向量
                eigenvectors = Xk1 * Phi_q;
                return { current_eigenvalues, eigenvectors };
            }

            // 更新迭代矩阵
            Xk = Xk1 * Phi;
            prev_eigenvalues = current_eigenvalues;
        }

        throw runtime_error("子空间迭代未收敛");
    }
   
    Matrix getColumn(int col) const {
        Matrix vec(rows, 1);
        for (int i = 0; i < rows; ++i) {
            vec(i, 0) = data[i][col];
        }
        return vec;
    }

    void setColumn(int col, const Matrix& vec) {
        if (vec.rows != rows || vec.cols != 1) {
            throw invalid_argument("输入必须是列向量");
        }
        for (int i = 0; i < rows; ++i) {
            data[i][col] = vec(i, 0);
        }
    }
    // 矩阵信息获取
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    //pair<Matrix, Matrix> StandardEigenSolver(const Matrix& A);
   
};




int main() {
    Pipe p1;
    
    //将读入数据放进p1
    // Matrix StiffnessMatrix(2, 2);
    // Matrix MassMatrix(2, 2);

    vector<int> id;       // 存储首次出现的数字
    int id_start=0, id_end=0;
    unordered_map<int, int> number; // 记录数字首次出现的位置（1-based）<编号，序号>

    string filename = "Default.txt";
    std::ifstream file(filename);
    json j;
    file >> j; 

    const auto& pipes = j["pipe_system"]["pipes"];

    
    int Num_Pipes = 0;
    for (const auto& pipe : pipes){const auto& elements = pipe["elements"];
        for (const auto& elem : elements) {Num_Pipes++;}}
    std::cout << "Num_Pipes: " << Num_Pipes << std::endl;  
    
    Matrix StiffnessMatrix(Num_Pipes*12, Num_Pipes*12);
    Matrix MassMatrix(Num_Pipes*12, Num_Pipes*12);
    // int id = 0;
    int Count_pipes2 = 0;

    for (const auto& pipe : pipes) {
        const auto& elements = pipe["elements"];
        for (const auto& elem : elements) {
            Count_pipes2++;
            p1.StartNode = elem.value("from_node", -1);
            p1.EndNode = elem.value("to_node", -1);
            p1.DeltaX = elem.value("delta_x", 0.0)/1000;
            p1.DeltaY = elem.value("delta_y", 0.0)/1000;
            p1.DeltaZ = elem.value("delta_z", 0.0)/1000;
            p1.Diameter = elem.value("diameter", 0.0)/1000;
            p1.WallThickness= elem.value("wall_thickness", 0.0)/1000;
            //double insulation_thickness = elem.value("insulation_thickness", 0.0);
            p1 .ElasticModulus= elem.value("elastic_modulus_cold", 0.0);
            if(p1 .ElasticModulus==0) p1 .ElasticModulus=200000000000;
            p1.PoissionRatio= elem.value("poisson_ratio", 0.0);
            if(p1.PoissionRatio==0) p1.PoissionRatio=0.3;

            // 读取 fluid_density 数组中的第一个值
            p1.Density= 0.0;
            if (elem.contains("fluid_density") && elem["fluid_density"].is_array() && !elem["fluid_density"].empty()) {
                p1.Density = elem["fluid_density"][0].get<double>();
            }



        p1.ShearModulus = p1.ElasticModulus / 2 / (1 + p1.PoissionRatio);
        p1.Length = pow((pow(p1.DeltaX, 2) + pow(p1.DeltaY, 2) + pow(p1.DeltaZ, 2)), 0.5);
        p1.Aera = Pi / 4 * (pow(p1.Diameter, 2) - pow((p1.Diameter - 2 * p1.WallThickness), 2));
        p1.Mass = p1.Aera * p1.Density * p1.Length;
        p1.I_x = Pi / 32 * (pow(p1.Diameter, 4) - pow((p1.Diameter - 2 * p1.WallThickness), 4)) * p1.Density * p1.Length;
        p1.J_x = Pi / 32 * (pow(p1.Diameter, 4) - pow((p1.Diameter - 2 * p1.WallThickness), 4));
        p1.J_y = p1.J_x / 2;


        //集中质量矩阵
        for (int i = 0; i <= 8; i++) {
            if (i <= 2) M1[i][i] = p1.Mass / 2;
            if (i >= 6) M1[i][i] = p1.Mass / 2;
        }
        /*//一致质量矩阵
        M1[0][0] = p1.Mass / 3; M1[6][0] = p1.Mass / 6;
        M1[1][1] = p1.Mass * 13 / 35; M1[5][1] = p1.Mass*p1.Length * 11 / 210; M1[7][1] = p1.Mass * 9 / 70; M1[11][1] = -p1.Mass * p1.Length * 13 / 420;
        M1[2][2] = p1.Mass * 13 / 35; M1[4][2] = -p1.Mass * p1.Length * 11 / 210; M1[8][2] = p1.Mass * 9 / 70; M1[10][2] = p1.Mass * p1.Length * 13 / 420;
        M1[3][3] = p1.I_x / 3; M1[9][3] = p1.I_x / 6;
        M1[4][4] = p1.Mass / 105 * pow(p1.Length, 2); M1[8][4] = -13 * p1.Mass / 420 * p1.Length; M1[10][4] = -p1.Mass / 140 * pow(p1.Length, 2);
        M1[5][5] = p1.Mass / 105 * pow(p1.Length, 2); M1[7][5] = 13 * p1.Mass / 420 * p1.Length; M1[11][5] = -p1.Mass / 140 * pow(p1.Length, 2);
        M1[6][6] = p1.Mass / 3;
        M1[7][7] = p1.Mass * 13 / 35; M1[11][7] = -p1.Mass * p1.Length * 11 / 210;
        M1[8][8] = p1.Mass * 13 / 35; M1[10][8] = -p1.Mass * p1.Length * 11 / 210;
        M1[9][9] = p1.I_x / 3;
        M1[10][10] = p1.Mass / 105 * pow(p1.Length, 2);
        M1[11][11] = p1.Mass / 105 * pow(p1.Length, 2);
        symmetry(M1,12);
        */

        //刚度矩阵
        K1[0][0] = p1.ElasticModulus * p1.Aera / p1.Length; K1[6][0] = -p1.ElasticModulus * p1.Aera / p1.Length;
        K1[1][1] = 12 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 3); K1[5][1] = 6 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 2); K1[7][1] = -12 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 3); K1[11][1] = 6 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 2);
        K1[2][2] = 12 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 3); K1[4][2] = -6 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 2); K1[8][2] = -12 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 3); K1[10][2] = 6 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 2);
        K1[3][3] = p1.J_x * p1.ShearModulus / p1.Length; K1[9][3] = -p1.J_x * p1.ShearModulus / p1.Length;
        K1[4][4] = 4 * p1.ElasticModulus * p1.J_y / p1.Length; K1[8][4] = 6 * p1.ElasticModulus * p1.J_y / p1.Length; K1[10][4] = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1[5][5] = 4 * p1.ElasticModulus * p1.J_y / p1.Length; K1[7][5] = -6 * p1.ElasticModulus * p1.J_y / p1.Length; K1[11][5] = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1[6][6] = p1.ElasticModulus * p1.Aera / p1.Length;
        K1[7][7] = 12 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 3); K1[11][7] = -6 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 2);
        K1[8][8] = 12 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 3); K1[10][8] = 6 * p1.ElasticModulus * p1.J_y / pow(p1.Length, 2);
        K1[9][9] = p1.J_x * p1.ShearModulus / p1.Length;
        K1[10][10] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1[11][11] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        symmetry(K1, 12);

        float h = pow(pow(p1.DeltaX, 2) + pow(p1.DeltaZ, 2), 0.5);
        if(h==0){
        c[0][0] = p1.DeltaX / p1.Length; c[0][1] = p1.DeltaY / p1.Length; c[0][2] = p1.DeltaZ / p1.Length;
        c[1][0] = 0; c[1][1] = h / p1.Length; c[1][2] = 1;
        c[2][0] = 1;c[2][1]=0; c[2][2] = 0;


        }
        else{
            
        c[0][0] = p1.DeltaX / p1.Length; c[0][1] = p1.DeltaY / p1.Length; c[0][2] = p1.DeltaZ / p1.Length;
        c[1][0] = -p1.DeltaX * p1.DeltaY / p1.Length / h; c[1][1] = h/p1.Length; c[1][2] = -p1.DeltaZ * p1.DeltaY / p1.Length / h;
        c[2][0] = -p1.DeltaZ / h;c[2][1]=0; c[2][2] = p1.DeltaX / h;
        }

        cout<< "c: " << endl;
        for (int i = 0; i < 3; i++) {  
            for (int j = 0; j < 3; j++) {
                cout << c[i][j] << " ";
            }
            cout << endl;
        }

        for (int i = 0; i <= 3; i++) {
            for (int j = 0; j <= 2; j++) {
                for (int k = 0; k <= 2; k++) {
                    C[4 * i + j][4 * i + k] = c[j][k];
                }
            }
        }
        for (int i = 0; i <= 11; i++) {
            for (int j = 0; j <= 11; j++) {
                CT[i][j] = C[j][i];
            }
        }

        if (number.find(p1.StartNode) == number.end()) {
            // 数字未出现过，记录位置并添加到vector
            number[p1.StartNode] = id.size() + 1; // 位置从1开始
            id.push_back(p1.StartNode);
        }
        // number[num]*6做矩阵组装
        id_start = number[p1.StartNode]-1;

        if (number.find(p1.EndNode) == number.end()) {
            // 数字未出现过，记录位置并添加到vector
            number[p1.EndNode] = id.size() + 1; // 位置从1开始
            id.push_back(p1.EndNode);
        }
        // number[num]*6做矩阵组装
        id_end = number[p1.EndNode]-1;

        //矩阵乘法：M=CT*M1*C,K=CT*K1*C，完成矩阵生成
        Multi(CT, M1, C);
        cout << "CT: " << endl;
        for (int i = 0; i <= 11; i++) {
            for (int j = 0; j <= 11; j++) {
                cout << CT[i][j] << " ";
            }
            cout << endl;
        }
        cout<<endl;
        cout << "M1: " << endl;
        for (int i = 0; i <= 11; i++) {
            for (int j = 0; j <= 11; j++) {
                cout << M1[i][j] << " ";
            }
            cout << endl;
        }
        cout<<endl;
        cout << "C: " << endl;
        for (int i = 0; i <= 11; i++) {
            for (int j = 0; j <= 11; j++) {
                cout << C[i][j] << " ";
            }
            cout << endl;
        }
        cout<<endl;

        
        for (int i = 0; i <= 5; i++) {
            for (int j = 0; j <= 5; j++) {
                MassMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult[i][j];
            }
        }

        for (int i = 6; i <= 11; i++) {
            for (int j = 0; j <= 5; j++) {
                MassMatrix(i + 6 * id_end-6, j + 6 * id_start) += MultiResult[i][j];
            }
        }

        for (int i = 0; i <= 5; i++) {
            for (int j = 6; j <= 11; j++) {
                MassMatrix(i + 6 * id_start, j + 6 * id_end-6) += MultiResult[i][j];
            }
        }

        for (int i = 6; i <= 11; i++) {
            for (int j = 6; j <= 11; j++) {
                MassMatrix(i + 6 * id_end-6, j + 6 * id_end-6) += MultiResult[i][j];
            }
        }

        Multi(CT, K1, C);

        for (int i = 0; i <= 5; i++) {
            for (int j = 0; j <= 5; j++) {
                StiffnessMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult[i][j];
            }
        }

        for (int i = 6; i <= 11; i++) {
            for (int j = 0; j <= 5; j++) {
                StiffnessMatrix(i + 6 * id_end-6, j + 6 * id_start) += MultiResult[i][j];
            }
        }

        for (int i = 0; i <= 5; i++) {
            for (int j = 6; j <= 11; j++) {
                StiffnessMatrix(i + 6 * id_start, j + 6 * id_end-6) += MultiResult[i][j];
            }
        }

        for (int i = 6; i <= 11; i++) {
            for (int j = 6; j <= 11; j++) {
                StiffnessMatrix(i + 6 * id_end-6, j + 6 * id_end-6) += MultiResult[i][j];
            }
        }


        //矩阵组装后可丢弃：
        memset(K1, 0, sizeof(K1)); memset(M1, 0, sizeof(M1)); memset(K, 0, sizeof(K)); memset(M, 0, sizeof(M)); memset(C, 0, sizeof(C)); memset(c, 0, sizeof(c)); memset(CT, 0, sizeof(CT));
        memset(&p1, 0, sizeof(p1));
        }
    }

    StiffnessMatrix.rows=id.size()*6;StiffnessMatrix.cols=id.size()*6;
    MassMatrix.rows=id.size()*6;MassMatrix.cols=id.size()*6;

    Matrix temp1(1,1);
    auto result = temp1.subspaceIteration(StiffnessMatrix,MassMatrix,8);
    auto eigenvalues = result.first;
    auto eigenvectors = result.second;
    // 遍历输出 eigenvalues 中的每个元素
    cout << "特征值: ";
    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        cout << eigenvalues[i];
        if (i < eigenvalues.size() - 1) {
            cout << ", ";
        }
    }
    cout << endl;
    // cout << eigenvectors;

	return 0;
}
