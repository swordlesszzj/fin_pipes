#include "Matrix.h"
#include "lapacke.h"
//#include <lapacke.h>
// #include <random> 
bool Matrix::isSquare() const
{
    return rows == cols;
}

Matrix::Matrix(int r, int c) : rows(r), cols(c)
{
    data.resize(rows, vector<double>(cols, 0.0));
}

Matrix::Matrix(const vector<vector<double>> &vec)
{
    if (vec.empty() || vec[0].empty())
    {
        throw invalid_argument("无效矩阵数据");
    }
    rows = vec.size();
    cols = vec[0].size();
    for (const auto &row : vec)
    {
        if (row.size() != cols)
        {
            throw invalid_argument("矩阵行长度不一致");
        }
    }
    data = vec;
}

Matrix Matrix::operator+(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
    {
        throw invalid_argument("矩阵维度不匹配");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result(i, j) = data[i][j] + other(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
    {
        throw invalid_argument("矩阵维度不匹配");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result(i, j) = data[i][j] - other(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix &other) const
{
    if (cols != other.rows)
    {
        throw invalid_argument("矩阵维度不匹配");
    }
    Matrix result(rows, other.cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < other.cols; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < cols; ++k)
            {
                sum += data[i][k] * other(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const
{
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result(i, j) = data[i][j] * scalar;
        }
    }
    return result;
}

Matrix operator*(double scalar, const Matrix &m)
{
    return m * scalar;
}

Matrix Matrix::operator/(double scalar) const
{
    if (scalar == 0)
    {
        throw invalid_argument("除数不能为零");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result(i, j) = data[i][j] / scalar;
        }
    }
    return result;
}

double &Matrix::operator()(int i, int j)
{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
    {
        throw out_of_range("索引越界");
    }
    return data[i][j];
}

const double &Matrix::operator()(int i, int j) const
{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
    {
        throw out_of_range("索引越界");
    }
    return data[i][j];
}

ostream &operator<<(ostream &os, const Matrix &m)
{
    for (int i = 0; i < m.rows; ++i)
    {
        for (int j = 0; j < m.cols; ++j)
        {
            os << m.data[i][j];
            if (j < m.cols - 1)
                os << "\t";
        }
        os << endl;
    }
    return os;
}

Matrix Matrix::transpose() const
{
    vector<vector<double>> transposed(cols, vector<double>(rows));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            transposed[j][i] = data[i][j];
        }
    }
    Matrix temp(transposed);
    return temp;
}

Matrix Matrix::Solve_Matrix(const Matrix &A, const Matrix &B)
{
    if (!A.isSquare())
    {
        throw invalid_argument("系数矩阵A必须为方阵");
    }
    if (A.getRows() != B.getRows())
    {
        throw invalid_argument("矩阵A和B行数必须匹配");
    }

    int n = A.getRows();
    int m = B.getCols();
    Matrix x(n, m);

    vector<vector<double>> Ab(n, vector<double>(n + m));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            Ab[i][j] = A(i, j);
        }
        for (int j = 0; j < m; ++j)
        {
            Ab[i][n + j] = B(i, j);
        }
    }

    const double EPS = 1e-40;
    for (int k1 = 0; k1 < n; ++k1)
    {
        int maxRow = k1;
        for (int i = k1; i < n; ++i)
        {
            if (fabs(Ab[i][k1]) > fabs(Ab[maxRow][k1]))
            {
                maxRow = i;
            }
        }

        //cout <<k1<<"   "<< fabs(Ab[maxRow][k1]) << endl;
        if (fabs(Ab[maxRow][k1]) < EPS)
        {
            
            throw domain_error("矩阵奇异，无法求解");
        }
        swap(Ab[k1], Ab[maxRow]);
        double pivot = Ab[k1][k1];
        for (int j = k1; j < n + m; ++j)
        {
            Ab[k1][j] /= pivot;
        }
        for (int i = 0; i < n; ++i)
        {
            if (i != k1 && fabs(Ab[i][k1]) > EPS)
            {
                double factor = Ab[i][k1];
                for (int j = k1; j < n + m; ++j)
                {
                    Ab[i][j] -= factor * Ab[k1][j];
                }
            }
        }
    }

    for (int col = 0; col < m; ++col)
    {
        for (int i = 0; i < n; ++i)
        {
            x(i, col) = Ab[i][n + col];
            for (int j = 0; j < i; ++j)
            {
                x(i, col) -= Ab[i][j] * x(j, col);
            }
        }
    }

    return x;
}

pair<vector<double>, Matrix> Matrix::lanczosIteration(const Matrix &K, const Matrix &M, int q, double tol)
{
    int n = K.rows;
    int r = min(2 * q, q + 8);

    Matrix x(n, 1); // 初始全1向量
    for (int i = 0; i < n; ++i)
    {
        x(i, 0) = 1.0;
    }

    Matrix temp = x.transpose() * M * x;
    double beta = sqrt(temp(0, 0));

    x = x / beta;

    Matrix X(n, r);    // Lanczos 向量矩阵
    X.setColumn(0, x); // 第一列

    vector<double> alpha(r, 0.0);
    vector<double> beta_vals(r, 0.0);
    beta_vals[0] = beta;
    Matrix x_hat(1,r);
    for (int i = 1; i < r; ++i)
    {
        Matrix b = M * X.getColumn(i - 1);
        x_hat = Solve_Matrix(K, b);

        temp = x_hat.transpose() * M * X.getColumn(i - 1);
        alpha[i - 1] = temp(0, 0);

        if (i == 1)
        {
            x_hat = x_hat - alpha[i - 1] * X.getColumn(i - 1);
        }
        else
        {
            x_hat = x_hat - alpha[i - 1] * X.getColumn(i - 1) - beta_vals[i - 1] * X.getColumn(i - 2);
        }

        temp = x_hat.transpose() * M * x_hat;
        if (temp(0, 0) < 0)
        {
            //throw invalid_argument("Lanczos迭代中出现负数");
            temp(0, 0)=-temp(0, 0);
            cout<<"Lanczos迭代中出现负数";
        }
        
        beta_vals[i] = sqrt(temp(0, 0));

        if (beta_vals[i] < tol)
        {
            r = i;
            throw invalid_argument("beta小于1e-40");
            break;
        }

        x_hat = x_hat / beta_vals[i];
        X.setColumn(i, x_hat);
    }
    x_hat = Solve_Matrix(K, M*x_hat);
    temp = x_hat.transpose() * M * x_hat;
    alpha[r - 1] = temp(0, 0);

    SymmetricTridiagonalQR T(r);
    for (int i = 0; i < r; ++i)
    {
        T.setDiagonal(i, alpha[i]);
        if (i < r - 1)
        {
            T.setOffDiagonal(i, beta_vals[i + 1]);
        }
    }

    auto result = T.computeEigenvaluesAndVectors();
    auto eigenvalues = result.first;
    auto Z = result.second;
    cout<<"Lanczos三对角矩阵特征值：";
    for(const auto &val: eigenvalues){
        cout<<val<<" ";
    }
    cout<<endl;
    for (int i=0;i<r;i++){
     cout<<"alpha"<<i+1<<"="<<alpha[i]<<"  ";
     if(i<r-1) cout<<"beta"<<i+1<<"="<<beta_vals[i+1];
     cout<<endl;
    }
    // for (const auto &val : eigenvalues
    if(eigenvalues.size()!=r) cout<<"三对角矩阵特征值少";
    Matrix Phi = X * Z;

    vector<double> omega(r);
    for (int i = 0; i < r; ++i)
    {
        omega[i] =eigenvalues[i]; //1.0 / sqrt(eigenvalues[i]);
    }

    //vector<double> omega_q(q);
    // Matrix Phi_q(n, q);

    // vector<size_t> indices(r);
    // for (size_t i = 0; i < r; ++i)
    //     indices[i] = i;
    // sort(indices.begin(), indices.end(), [&](size_t a, size_t b)
    //      { return omega[a] < omega[b]; });

    // for (int i = 0; i < q; ++i)
    // {
    //     omega_q[i] = omega[indices[i]];
    //     Phi_q.setColumn(i, Phi.getColumn(indices[i]));
    // }

    //cout<<"里兹法频率：";
    //for(int i=0;i<q;i++) cout<<omega_q[i]<<" ";
    //cout<<endl;

    return make_pair(omega, Phi);
}

Matrix Matrix::choleskyDecomposition(const Matrix &M)
{
    int n = M.rows;
    Matrix L(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < j; ++k)
            {
                sum += L(i, k) * L(j, k);
            }
            if (i == j)
            {
                double diag = M(i, i) - sum;
                double EPS=1e-8;
                if(diag<0) {cout<<"非正定";diag=fabs(M(i,i));}//fabs(diag/M(i, i)<EPS)||
                if(fabs(diag/M(i, i))<EPS){cout<<"半正定";diag=fabs(M(i,i));}
                if (diag < 0)
                {
                    std::cerr << "Cholesky分解失败：矩阵不是正定的，M(" << i << "," << i << ") - sum = " << diag << std::endl;
                    throw std::runtime_error("Cholesky分解失败：矩阵不是正定的");
                }
                L(i, j) = sqrt(diag);
            }
            else
            {
            //    if (fabs(L(j, j)) < 1e-58)
            //         L(i, j) = 0.0; // 对应主元为0，下三角也为0
            //     else
                    L(i, j) = (M(i, j) - sum) / L(j, j);
            }
        }
    }
    return L;
}


pair<vector<double>, Matrix> Matrix::subspaceIteration(const Matrix &K, const Matrix &M, int q, double tol, int max_iter)
{
    int n = K.rows;
    int r = min(2 * q, q + 8);
    //Matrix ML=choleskyDecomposition(M);
    //Matrix MM=ML*ML.transpose();
    // double max=0,max1=0;
    // for(int i=0;i<M.rows;i++){
    //     for(int j=0;j<M.cols;j++){
    //         if(fabs(M(i,j))>max) max=fabs(M(i,j));
    //         if(fabs(MM(i,j)-M(i,j))>max1) max1=fabs(MM(i,j)-M(i,j));
    //     }
    // }
    // cout<<max<<"  "<<max1;
    // Matrix X0(n, r);
    // std::mt19937 gen(42); // 固定种子，便于复现
    // std::normal_distribution<> dist(0.0, 1.0);
    // for (int i = 0; i < n; ++i)
    //     for (int j = 0; j < r; ++j)
    //         X0(i, j) = dist(gen);

    // // Gram-Schmidt 正交化
    // for (int j = 0; j < r; ++j) {
    //     // 对前j-1列做投影
    //     for (int k = 0; k < j; ++k) {
    //         double dot = 0.0, norm = 0.0;
    //         for (int i = 0; i < n; ++i) {
    //             dot += X0(i, j) * X0(i, k);
    //             norm += X0(i, k) * X0(i, k);
    //         }
    //         if (norm > 1e-12) {
    //             for (int i = 0; i < n; ++i)
    //                 X0(i, j) -= dot / norm * X0(i, k);
    //         }
    //     }
    //     // 归一化
    //     double norm = 0.0;
    //     for (int i = 0; i < n; ++i)
    //         norm += X0(i, j) * X0(i, j);
    //     norm = sqrt(norm);
    //     if (norm > 1e-12) {
    //         for (int i = 0; i < n; ++i)
    //             X0(i, j) /= norm;
    //     }
    // }
    
    Matrix X0(n, r);
    for (int i = 0; i < n; ++i)
    {
        X0(i, 0) = 1.0;
    }
    for (int i = 1; i < r; ++i)
    {
        if (i - 1 < n)
            X0(i - 1, i) = 1.0;
    }

    Matrix Xk = X0;
    vector<double> prev_eigenvalues(q, 0.0);
    Matrix eigenvectors(n, q);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        cout << "迭代次数: " << iter + 1 << endl;
        //    cout<<M<<endl;
        //    cout<<Xk<<endl;
        Matrix Y = M * Xk;

        Matrix Xk1 = Solve_Matrix(K, Y);

        Matrix K_star = Xk1.transpose() * Y;//Xk1.transpose() *K*Xk1;

        //Matrix Y1 = ML* Xk1;

        Matrix M_star = Xk1.transpose() *M*Xk1;

        Matrix L =choleskyDecomposition(M_star);
        //cout<<K_star<<endl<<M_star;

        // auto result = lanczosIteration(K_star, M_star, q);
        // auto eigenvalues = result.first;
        // auto Z = result.second;
        int n = K_star.getRows();
        std::vector<double> a(n * n);
        for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        a[i * n + j] = K_star(i, j); // 行主序

        std::vector<double> w(n); // 存储特征值
        int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, a.data(), n, w.data());
        if (info != 0) {
            throw std::runtime_error("LAPACKE_dsyev 计算特征值失败");
        }

        // a.data() 现在存储特征向量，每列一个
        Matrix Z(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                Z(i, j) = a[j * n + i]; // 注意LAPACK输出为列主序

        std::vector<double> eigenvalues(w.begin(), w.end());


        Matrix Phi = Solve_Matrix(L.transpose(), Z);

        vector<double> current_eigenvalues(r);
        Matrix Phi_q(r, r);
        for (int i = 0; i < r; ++i)
        {
            current_eigenvalues[i] = eigenvalues[i];
            for (int j = 0; j < r; ++j)
            {
                Phi_q(j, i) = Phi(j, i);
            }
        }
        cout<<"当前特征值：";
        for (const auto &val : current_eigenvalues)
        {
            cout << val << " ";
        }
        cout << endl;
        bool converged = true;
        for (int i = 0; i < r; ++i)
        {
            if (abs(current_eigenvalues[i] - prev_eigenvalues[i]) > tol * abs(current_eigenvalues[i]))
            {
                converged = false;
                break;
            }
        }

        if (converged)
        {
            eigenvectors = Xk1 * Phi_q;
            return {current_eigenvalues, eigenvectors};
        }
        double length;
        

        Xk = Xk1 * Phi;

        prev_eigenvalues = current_eigenvalues;
    }

    throw runtime_error("子空间迭代未收敛");
}

Matrix Matrix::getColumn(int col) const
{
    Matrix vec(rows, 1);
    for (int i = 0; i < rows; ++i)
    {
        vec(i, 0) = data[i][col];
    }
    return vec;
}

void Matrix::setColumn(int col, const Matrix &vec)
{
    if (vec.rows != rows || vec.cols != 1)
    {
        throw invalid_argument("输入必须是列向量");
    }
    for (int i = 0; i < rows; ++i)
    {
        data[i][col] = vec(i, 0);
    }
}

int Matrix::getRows() const
{
    return rows;
}

int Matrix::getCols() const
{
    return cols;
}
