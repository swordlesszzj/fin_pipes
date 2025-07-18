#include "Matrix.h"

bool Matrix::isSquare() const {
    return rows == cols;
}

Matrix::Matrix(int r, int c) : rows(r), cols(c) {
    data.resize(rows, vector<double>(cols, 0.0));
}

Matrix::Matrix(const vector<vector<double>>& vec) {
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

Matrix Matrix::operator+(const Matrix& other) const {
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

Matrix Matrix::operator-(const Matrix& other) const {
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

Matrix Matrix::operator*(const Matrix& other) const {
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

Matrix Matrix::operator*(double scalar) const {
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = data[i][j] * scalar;
        }
    }
    return result;
}

Matrix operator*(double scalar, const Matrix& m) {
    return m * scalar;
}

Matrix Matrix::operator/(double scalar) const {
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

double& Matrix::operator()(int i, int j) {
    if (i < 0 || i >= rows || j < 0 || j >= cols) {
        throw out_of_range("索引越界");
    }
    return data[i][j];
}

const double& Matrix::operator()(int i, int j) const {
    if (i < 0 || i >= rows || j < 0 || j >= cols) {
        throw out_of_range("索引越界");
    }
    return data[i][j];
}

ostream& operator<<(ostream& os, const Matrix& m) {
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j) {
            os << m.data[i][j];
            if (j < m.cols - 1) os << "\t";
        }
        os << endl;
    }
    return os;
}

Matrix Matrix::transpose() const {
    vector<vector<double>> transposed(cols, vector<double>(rows));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = data[i][j];
        }
    }
    Matrix temp(transposed);
    return temp;
}

Matrix Matrix::Solve_Matrix(const Matrix& A, const Matrix& B) {
    if (!A.isSquare()) {
        throw invalid_argument("系数矩阵A必须为方阵");
    }
    if (A.getRows() != B.getRows()) {
        throw invalid_argument("矩阵A和B行数必须匹配");
    }

    int n = A.getRows();
    int m = B.getCols();
    Matrix x(n, m);

    vector<vector<double>> Ab(n, vector<double>(n + m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Ab[i][j] = A(i, j);
        }
        for (int j = 0; j < m; ++j) {
            Ab[i][n + j] = B(i, j);
        }
    }

    const double EPS = 1e-8;
    for (int k = 0; k < n; ++k) {
        int maxRow = k;
        for (int i = k; i < n; ++i) {
            if (abs(Ab[i][k]) > abs(Ab[maxRow][k])) { maxRow = i; }
        }
        if (abs(Ab[maxRow][k]) < EPS) {
            throw domain_error("矩阵奇异，无法求解");
        }
        swap(Ab[k], Ab[maxRow]);
        double pivot = Ab[k][k];
        for (int j = k; j < n + m; ++j) {
            Ab[k][j] /= pivot;
        }
        for (int i = 0; i < n; ++i) {
            if (i != k && abs(Ab[i][k]) > EPS) {
                double factor = Ab[i][k];
                for (int j = k; j < n + m; ++j) {
                    Ab[i][j] -= factor * Ab[k][j];
                }
            }
        }
    }

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

pair<vector<double>, Matrix> Matrix::lanczosIteration(const Matrix& K,const Matrix& M,int q,double tol) {
    int n = K.rows; 
    int r = min(2 * q, q + 8);  

    Matrix x(n, 1);  // 初始全1向量
    for (int i = 0; i < n; ++i) {
        x(i, 0) = 1.0;  
    }

    Matrix temp= x.transpose() * M * x; 
    double beta = sqrt(temp(0,0));
   
    x = x / beta;

    Matrix X(n, r);  // Lanczos 向量矩阵
    X.setColumn(0, x);  // 第一列

    vector<double> alpha(r, 0.0); 
    vector<double> beta_vals(r , 0.0);
    beta_vals[0] = beta;

    for (int i = 1; i < r; ++i) {
        Matrix b = M * X.getColumn(i - 1);
        Matrix x_hat = Solve_Matrix(K, b);

        temp= x_hat.transpose() * M * X.getColumn(i - 1);
        alpha[i - 1] = temp(0, 0);
        
        if (i == 1) {
            x_hat = x_hat - alpha[i - 1] * X.getColumn(i - 1);
        }
        else {
            x_hat = x_hat - alpha[i - 1] * X.getColumn(i - 1)
                - beta_vals[i - 1] * X.getColumn(i - 2);
        }

        temp = x_hat.transpose() * M * x_hat;
        beta_vals[i] = sqrt(temp(0,0));

        if (beta_vals[i] < tol) {
            r = i;
            break;
        }
        
        x_hat = x_hat / beta_vals[i];
        X.setColumn(i, x_hat); 
    }

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

    Matrix Phi = X * Z;

    vector<double> omega(r);
    for (int i = 0; i < r; ++i) {
        omega[i] = 1.0 / sqrt(eigenvalues[i]);
    }

    vector<double> omega_q(q);
    Matrix Phi_q(n, q);

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

Matrix Matrix::choleskyDecomposition(const Matrix& M) {
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

pair<vector<double>, Matrix> Matrix::subspaceIteration(const Matrix& K, const Matrix& M, int q, double tol, int max_iter) {
    int n = K.rows; 
    int r = min(2 * q, q + 8); 

    Matrix X0(n, r);
    for (int i = 0; i < n; ++i) {
        X0(i, 0) = 1.0;
    }
    for (int i = 1; i < r; ++i) {
        if (i - 1 < n) X0(i - 1, i) = 1.0;
    }

    Matrix Xk = X0;  
    vector<double> prev_eigenvalues(q, 0.0);  
    Matrix eigenvectors(n, q);  
    
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

        auto result=lanczosIteration(K_star, M_star,r);
        auto eigenvalues = result.first;
        auto Z = result.second;

        Matrix Phi = Solve_Matrix(L.transpose(), Z);

        vector<double> current_eigenvalues(q);
        Matrix Phi_q(r, q);
        for (int i = 0; i < q; ++i) {
            current_eigenvalues[i] = eigenvalues[i];
            for (int j = 0; j < r; ++j) {
                Phi_q(j, i) = Phi(j, i);
            }
        }

        bool converged = true;
        for (int i = 0; i < q; ++i) {
            if (abs(current_eigenvalues[i] - prev_eigenvalues[i]) > tol * abs(current_eigenvalues[i])) {
                converged = false;
                break;
            }
        }

        if (converged) {
            eigenvectors = Xk1 * Phi_q;
            return { current_eigenvalues, eigenvectors };
        }

        Xk = Xk1 * Phi;
        prev_eigenvalues = current_eigenvalues;
    }

    throw runtime_error("子空间迭代未收敛");
}

Matrix Matrix::getColumn(int col) const {
    Matrix vec(rows, 1);
    for (int i = 0; i < rows; ++i) {
        vec(i, 0) = data[i][col];
    }
    return vec;
}

void Matrix::setColumn(int col, const Matrix& vec) {
    if (vec.rows != rows || vec.cols != 1) {
        throw invalid_argument("输入必须是列向量");
    }
    for (int i = 0; i < rows; ++i) {
        data[i][col] = vec(i, 0);
    }
}

int Matrix::getRows() const {
    return rows;
}

int Matrix::getCols() const {
    return cols;
}
