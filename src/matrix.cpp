#include "../include/matrix.h"
#include "mslice.h"
using namespace std;

/* CONSTRUCTORS */
Matrix::Matrix(const string file, const string delim = ",", const bool has_header = true)
{
    ifstream in;
    in.open(file);
    string first_line;
    string line;

    //double value;
    if (in.is_open())
    {
        if (has_header)
            getline(in, first_line); // get the first line

        size_t row_cnt = 0;
        size_t *col_cnt = 0;
        size_t *ncol;

        while (getline(in, line))
        {
            if (line.empty())
            {
                continue;
                row_cnt -= 1;
            }
            split_str(line, delim);
            string token;
            istringstream st(line);
            size_t col_cnt = 0;
            ncol = &col_cnt;
            while (st >> token)
            {
                double value;
                try
                {
                    value = std::stod(token);
                    col_cnt += 1;
                }
                catch (const std::exception &e)
                {
                    //col_cnt -= 1;
                    continue;
                }
                mat_data.push_back(value);
            }
            row_cnt += 1;
        }

        nb_col = *ncol;
        nb_row = row_cnt;
    }
    else
    {
        std::cerr << "Couldn't open file" << std::endl;
    }
    in.close();
}

Matrix::Matrix(const size_t n, const size_t m, bind_d dist) : nb_row(n), nb_col(m)
{
    mat_data.resize(n * m);
    generate(mat_data.begin(), mat_data.end(), dist);
}

Matrix::Matrix(const size_t n, const size_t m, double fill = 0) : nb_row(n), nb_col(m)
{
    mat_data = vector_d(n * m, fill);
}

Matrix::Matrix(initializer_list<vector<double>> li)
{
    nb_row = li.size();
    auto b = li.begin();
    nb_col = b->size();
    mat_data.reserve(nb_row * nb_col);
    for (auto col : li)
    {
        if (nb_col != col.size())
            throw runtime_error("Matrix columns must have the same size !");
        for (double i : col)
            mat_data.emplace_back(i);
    }
}

/* COPY CONTROL */
Matrix::Matrix(const Matrix &cls)
{
    nb_row = cls.nb_row;
    nb_col = cls.nb_col;
    mat_data = cls.mat_data;
}

/* OPERATORS OVERLOADING */
Matrix &Matrix::operator=(const Matrix &rmat)
{
    nb_row = rmat.nb_row;
    nb_col = rmat.nb_col;
    mat_data = rmat.mat_data;
    return *this;
}

bool operator==(const Matrix &lmat, const Matrix &rmat)
{
    return (lmat.nb_row == rmat.nb_row) &&
           (lmat.nb_col == rmat.nb_col) &&
           (lmat.mat_data == rmat.mat_data);
}

bool operator!=(const Matrix &lmat, const Matrix &rmat)
{
    return !(lmat == rmat);
}

double &Matrix::operator()(const size_t n, const size_t m)
{
    return mat_data[idx(n, m, nb_col)];
}

double Matrix::operator()(const size_t n, const size_t m) const
{
    return mat_data[idx(n, m, nb_col)];
}

Mslice Matrix::operator()(const range row, const range col)
{
    Mslice M(row, col, this);
    return M;
}

double &Matrix::operator[](const size_t i) { return mat_data[i]; }

double Matrix::operator[](const size_t i) const { return mat_data[i]; }

Matrix &Matrix::operator+=(const Matrix &rmat)
{
    assert((nb_row == rmat.nb_row && nb_col == rmat.nb_col) && "Matrices with different dimensions");

    for (size_t i = 0; i != nb_row * nb_col; i++)
    {
        mat_data[i] += rmat[i];
    }
    return *this;
}

Matrix operator+(const Matrix &lmat, const Matrix &rmat)
{
    Matrix sum = lmat;
    sum += rmat;
    return sum;
}

Matrix operator+(const Matrix &lmat, const double val)
{
    Matrix sum = lmat;
    for (size_t n = 0; n != lmat.get_nb_rows() * lmat.get_nb_columns(); n++)
    {
        sum[n] += val;
    }
    return sum;
}

Matrix operator+(const double val, const Matrix &lmat)
{
    return lmat + val;
}

Matrix &Matrix::operator-=(const Matrix &rmat)
{
    assert((nb_row == rmat.nb_row && nb_col == rmat.nb_col) && "Matrices with different dimensions");
    for (size_t i = 0; i != nb_row * nb_col; i++)
    {
        mat_data[i] -= rmat[i];
    }
    return *this;
}

Matrix operator-(const Matrix &lmat, const Matrix &rmat)
{
    Matrix sum = lmat;
    sum -= rmat;
    return sum;
}

Matrix operator-(const Matrix &lmat, const double val)
{
    Matrix sum = lmat;
    for (size_t n = 0; n != lmat.get_nb_rows() * lmat.get_nb_columns(); n++)
    {
        sum[n] -= val;
    }
    return sum;
}

Matrix operator-(const double val, const Matrix &lmat)
{
    return lmat - val;
}

Matrix &Matrix::operator*=(const double val)
{
    for (size_t i = 0; i != nb_row * nb_col; i++)
    {
        mat_data[i] *= val;
    }
    return *this;
}

Matrix operator*(const Matrix &lmat, const double val)
{
    Matrix res = lmat;
    res *= val;
    return res;
}

Matrix operator*(const double val, const Matrix &lmat)
{
    return lmat * val;
}

Matrix operator*(const Matrix &lmat, const Matrix &rmat)
{
    size_t nrow = lmat.get_nb_rows();
    size_t mcol = lmat.get_nb_columns();
    size_t kcol = rmat.get_nb_columns();

    assert(mcol == rmat.get_nb_rows() && "Number of columns must match number of rows");
    Matrix res(nrow, kcol);
    for (size_t i = 0; i != nrow; i++)
    {
        for (size_t j = 0; j != mcol; j++)
        {
            for (size_t k = 0; k != kcol; k++)
            {
                res(i, k) += lmat(i, j) * rmat(j, k);
            }
        }
    }
    return res;
}

Matrix operator/(const Matrix &num_mat, const double val)
{
    return num_mat * (1.0 / val);
}

// divide a matrix by a row matrix
Matrix operator/(const Matrix &num_mat, const Matrix &den_mat)
{
    auto n = num_mat.dim().first;
    auto m = num_mat.dim().second;
    assert(den_mat.dim().second == m && den_mat.dim().first == 1);
    double d;
    Matrix R(n, m);
    for (size_t i = 0; i != m; i++)
    {
        for (size_t j = 0; j != n; j++)
        {
            d = den_mat(0, i);
            if (d == 0)
                d = 1;
            R(j, i) = num_mat(j, i) / d;
        }
    }
    return R;
}

/* FUNCTIONS for SLICING, SUBSCRIPTING, CONVERTING A MATRIX */

vector_d Matrix::column(size_t m, size_t starting_row = 0)
{
    assert(m < nb_col && "Out of range");
    vector_d vcol;
    vcol.reserve(nb_col - starting_row);
    for (size_t i = starting_row; i != nb_row; ++i)
    {
        vcol.emplace_back(mat_data[idx(i, m, nb_col)]); // start = m, stride = nb_col
    }
    return vcol;
}

vector_d Matrix::column(size_t m, size_t starting_row = 0) const
{
    assert(m < nb_col && "Out of range");
    vector_d vcol;
    vcol.reserve(nb_col - starting_row);
    for (size_t i = starting_row; i != nb_row; ++i)
    {
        vcol.emplace_back(mat_data[idx(i, m, nb_col)]); // start = m, stride = nb_col
    }
    return vcol;
}

vector_d Matrix::row(size_t n)
{
    assert(n < nb_row && "Out of range");
    vector_d vrow;
    vrow.reserve(nb_row);
    for (size_t i = 0; i != nb_col; ++i)
    {
        vrow.emplace_back(mat_data[idx(i, n * nb_col, 1)]);
    }

    return vrow;
}

vector_d Matrix::row(size_t n) const
{
    assert(n < nb_row && "Out of range");
    vector_d vrow;
    vrow.reserve(nb_row);
    for (size_t i = 0; i != nb_col; ++i)
    {
        vrow.emplace_back(mat_data[idx(i, n * nb_col, 1)]);
    }
    return vrow;
}

Matrix as_matrix(vector<double> x, size_t axe = 1)
{
    assert((axe == 1) || (axe == 0) && "for matrix row set axe=0, for matrix column set axe=1");
    Matrix M;
    M.mat_data = x;
    if (axe == 1)
    {
        M.nb_row = x.size();
        M.nb_col = 1;
    }
    else if (axe == 0)
    {
        M.nb_row = 1;
        M.nb_col = x.size();
    }
    return M;
}

vector_d as_vector(const Matrix &A)
{
    return A.mat_data;
}

Matrix Matrix::sub_matrix(const range row, const range col) const
{
  
    if ((row.end > nb_row) || (col.end > nb_col))
    {
        throw runtime_error("Out of range value");
    }

    Matrix M;
    M.nb_row = row.len;
    M.nb_col = col.len;
    M.mat_data.reserve(M.nb_row * M.nb_col);

    for (auto i = row.start; i != row.end; i++)
    {
        for (auto j = col.start; j != col.end; j++)
        {
            M.mat_data.emplace_back(mat_data[idx(i, j, nb_col)]);
        }
    }
    return M;
}

void Matrix::add_row(const vector<double> vect, const size_t pos)
{
    assert(vect.size() == nb_col && "The vector must has the same size as nb_col");
    vector<double>::iterator it;
    it = mat_data.begin() + pos * nb_col;
    mat_data.insert(it, vect.begin(), vect.end());
    nb_row += 1;
}

void Matrix::add_column(const vector<double> vect, const size_t pos)
{
    assert(vect.size() == nb_row && "The vector must has the same size as nb_rows");
    for (size_t i = 0; i != nb_row; i++)
    {
        vector<double>::iterator it;
        it = mat_data.begin() + (i * (nb_col + 1) + pos);
        mat_data.insert(it, vect[i]);
    }
    nb_col += 1;
}

void Matrix::remove_column(const size_t col)
{
    assert(col < nb_row && "Out of range");
    for (size_t r = 0; r != nb_row; r++)
    {
        mat_data.erase(mat_data.begin() + idx(r, col, nb_col - 1));
    }
    nb_col -= 1;
}

void Matrix::split_str(string &str, string delimiter)
{
    size_t pos;
    while ((pos = str.find(delimiter)) != std::string::npos)
    {
        str[pos] = ' ';
    }
}

bool desc(double a, double b)
{
    return (a > b);
}

//sort a flat matrix
void Matrix ::sort_matrix(string s = "descending")
{
    assert(s == "descending" || s == "ascending");
    if (s == "descending")
        sort(mat_data.begin(), mat_data.end(), desc);
    else
        sort(mat_data.begin(), mat_data.end());
}

//reorder columns
void Matrix::reorder_columns(vector<size_t> vec)
{
    vector<double> newdata;
    newdata.reserve(nb_col * nb_row);
    for (auto v : vec)
    {
        for (size_t i = 0; i != nb_row; ++i)
        {
            newdata.emplace_back(mat_data[idx(i, v, nb_col)]);
        }
    }
    mat_data = newdata;
    *this = this->T();
}

void Matrix::shape()
{
    std::cout << "(nb_rows, nb_cols) = "
              << "(" << nb_row << ", " << nb_col << ")" << std::endl;
}

void Matrix::reshape(const size_t n, const size_t m)
{
    assert(n * m == mat_data.size());
    nb_row = n;
    nb_col = m;
}

/* MATRIX TRANSFORMATIONS  */
Matrix Matrix::Id(const size_t n)
{
    Matrix id(n, n);

    for (size_t i = 0; i != n; i++)
    {
        id(i, i) = 1.0;
    }
    return id;
}

Matrix identity_mat(const size_t n)
{
    Matrix I;
    return I.Id(n);
}

Matrix diag_matrix(const Matrix &M)
{
    size_t nbcol = M.get_nb_columns();
    Matrix D(nbcol, nbcol);
    double val;
    for (size_t i = 0; i != nbcol; i++)
    {
        if (M.get_nb_rows() == 1)
        {
            D(i, i) = M(0, i);
        }
        else
        {
            D(i, i) = M(i, i);
        }
    }
    return D;
}

//return row matrix containing the diag
Matrix diag(const Matrix &M)
{
    size_t nbcol = M.get_nb_columns();
    vector_d d;
    d.reserve(nbcol);
    for (size_t i = 0; i != nbcol; i++)
    {
        d.emplace_back(M(i, i));
    }
    Matrix D = {d};
    return D;
};

Matrix Matrix ::T() const
{
    Matrix T;
    T.nb_row = nb_col;
    T.nb_col = nb_row;
    vector<double> newdata;
    newdata.reserve(nb_col * nb_row);
    for (size_t c = 0; c != nb_col; ++c)
    {
        for (size_t i = 0; i != nb_row; ++i)
        {
            newdata.emplace_back(mat_data[idx(i, c, nb_col)]);
        }
    }
    T.mat_data = newdata;
    return T;
}

Matrix transpose(const Matrix &M)
{
    return M.T();
}

double sgn(double val)
{
    return (val != 0) ? (double(0) < val) - (val < double(0)) : 1;
}

double norm(const vector<double> &vec)
{
    double n = 0;
    for (auto v : vec)
    {
        n += pow(v, 2);
    }
    return sqrt(n);
}

double norm(const Matrix &vec)
{
    return norm(as_vector(vec));
}

bool is_triangular(const Matrix &T)
{
    //double eps = std::numeric_limits<double>::epsilon();
    double eps = 0.0000001;

    for (size_t i = 1; i != T.dim().first; i++)
    {
        for (size_t j = 0; j != i; j++)
        {

            if (abs(T(i, j)) > eps)

                return false;
        }
    }
    return true;
}

/* FUNCTIONS & MEMBERS FOR STAT OPERATIONS */
Matrix elm2pow_n(const Matrix &M, double n)
{
    Matrix R;
    R.nb_row = M.nb_row;
    R.nb_col = M.nb_col;
    R.mat_data.resize(R.nb_elm());
    for (size_t i = 0; i != M.nb_elm(); i++)
    {
        R.mat_data[i] = pow(M.mat_data[i], n);
    }
    return R;
};

double sum(const Matrix &M)
{
    return M.sum();
}

Matrix tile(Matrix A, const size_t n)
{
    auto r = as_vector(A);
    for (size_t i = 1; i != n; i++)
        A.add_row(r, i);
    return A;
}

double avg(vector<double> v) { return accumulate(v.begin(), v.end(), 0.0) / v.size(); }

Matrix mean(const Matrix &A, const string dimension = "col")
{
    assert(dimension == "col" || dimension == "row" || dimension == "flat");
    Matrix means;
    vector_d m;
    if (dimension == "col")
    {
        for (size_t i = 0; i != A.dim().second; i++)
        {
            // cout<< avg(A.column(i));
            m.push_back(avg(A.column(i)));
        }
    }
    else if (dimension == "row")
    {
        for (size_t i = 0; i != A.dim().first; i++)
            m.push_back(avg(A.row(i)));
    }
    else if (dimension == "flat")
    {
        m.push_back(A.avg());
    }
    means = {m};
    return means;
}

Matrix stdev(const Matrix &M, const string dimension = "col")
{
    assert((dimension == "col") || (dimension == "row") && "dimension must be a string : col or row ");
    size_t N = M.nb_elm();
    size_t m = M.dim().second;
    size_t n = M.dim().first;
    Matrix R(1, m);
    Matrix B = mean(M, dimension);
    for (size_t j = 0; j != m; j++)
    {
        for (size_t i = 0; i != n; i++)
        {
            R(0, j) += pow(M(i, j) - B(0, j), 2);
        }
        if (dimension == "col")
            R(0, j) = sqrt(R(0, j) / n);
        else
            R(0, j) = sqrt(R(0, j) / m);
    }
    return R;
}

Matrix scale(const Matrix &M, const string dimension = "col")
{
    assert((dimension == "col") || (dimension == "row") && "dimension must be a string : col or row ");
    size_t m = M.dim().second;
    size_t n = M.dim().first;
    Matrix R(n, m);
    Matrix B = mean(M, dimension);
    Matrix S = stdev(M, dimension);
    for (size_t j = 0; j != m; j++)
    {
        auto s = S(0, j);
        if (s == 0)
            s = 1;
        for (size_t i = 0; i != n; i++)
        {
            R(i, j) = (M(i, j) - B(0, j)) / s;
        }
    }
    return R;
}

/* double round_val(double value, double nb_decimals = 2)
{
    auto x = pow(10, nb_decimals);
    return floor((value * x) + .5) / x;
} */

double uniform_dist(double lower, double upper)
{
    random_device rd;
    mt19937 mt(rd());
    std::uniform_real_distribution<> distribution(lower, upper);
    return distribution(mt);
}

double normal_dist(double mean, double stddev)
{
    random_device rd;
    mt19937 mt(rd());
    std::normal_distribution<> distribution(mean, stddev);
    return distribution(mt);
}

_Bind<double (*(double, double))(double, double)> unif(double lower, double upper)
{
    return bind(uniform_dist, lower, upper);
}

_Bind<double (*(double, double))(double, double)> rnorm(double mean, double stddev)
{
    return bind(normal_dist, mean, stddev);
}

/* ALGORITHMS */

/**
 * @brief Householder transformation of a matrix column
 * 
 * @param x is a column of a matrix
 * @return Matrix 
 */
Matrix householder(Matrix x)
{
    auto normx = norm(x);
    if (normx == 0)
    {
        x[0] = 1;
    }
    size_t dim = x.get_nb_rows();
    Matrix e(dim, 1);
    auto alpha = sgn(x[0]) * normx;
    e[0] = alpha;
    auto v = x + e;
    auto vT = transpose(v);
    auto vvT = v * vT;
    auto vTv = vT * v;
    auto H = identity_mat(dim) - ((2.0 * vvT) / vTv[0]);
    return H;
}

/**
 * @brief Transform a matrix to hessenberg form 
 * 
 * @param A The matrix to be transformed
 * @return Matrix type
 */
Matrix hessenberg(const Matrix &A)
{
    auto dim = A.dim().first;
    auto HAH = A;
    for (size_t i = 0; i != dim - 1; i++)
    {
        auto col = as_matrix(HAH.column(i, i + 1));
        auto hous = householder(col);
        auto H = identity_mat(dim);
        H(range(i + 1, dim - 1), range(i + 1, dim - 1)) = hous;
        HAH = H * HAH * H;
    }
    return HAH;
}

/**
 * @brief QR factorization of a matrix, A = QR
 * 
 * @param A The matrix to be factorized
 * @param Q The orthogonal matrix
 * @param R The triangular matrix
 */
void QR_factorization(Matrix A, Matrix *Q, Matrix *R)
{
    size_t dim;
    const size_t nrow = A.dim().first;
    const size_t ncol = A.dim().second;
    (ncol < nrow) ? dim = ncol : dim = nrow;

    *Q = identity_mat(nrow);
    *R = A;
    Matrix II;
    for (size_t i = 0; i != dim; i++)
    {
        auto col = as_matrix(A.column(0));
        auto hous = householder(col);
        if (i == 0)
            II = hous;
        else
        {
            II = identity_mat(nrow);
            II(range(i, nrow - 1), range(i, nrow - 1)) = hous;
        }
        *R = II * (*R);
        *Q = (*Q) * II;
        auto HA = hous * A;
        A = HA.sub_matrix(range(1, HA.get_nb_rows() - 1), range(1, HA.get_nb_columns() - 1));
    }
}

double wilkinson_shift(const double a, const double b, const double c)
{
    double eps = std::numeric_limits<double>::epsilon();
    double s = (a - c) / 2;
    return (c - sgn(s) * pow(b, 2)) / (abs(s) + sqrt(pow(s, 2) + pow(b, 2)));
}

/**
 * @brief QR algorithm with the wilkinson shift 
 * 
 * @param M The data
 * @param vals Store the matrix containing the eigenvalues
 * @param vects Store the matrix containing the eigenvectors
 */
void QR_algorithm_w_shift(const Matrix &M, Matrix *vals, Matrix *vects)
{
    auto dim = M.get_nb_columns();
    *vects = identity_mat(dim);
    Matrix I = identity_mat(dim);
    *vals = hessenberg(M);
    *vals = M;
    auto nb_iter = 0;
    while (!is_triangular(*vals))
    {
        // wilkinson shift
        double a, b, c, u;
        auto m = vects->get_nb_columns() - 1;
        a = (*vals)(m - 1, m - 1);
        b = (*vals)(m - 1, m);
        c = (*vals)(m, m);
        u = wilkinson_shift(a, b, c);

        // q, r, matrices
        Matrix q, r;
        QR_factorization(*vals - (u * I), &q, &r);
        *vals = r * q + (u * I);
        *vects = *vects * q;
        nb_iter++;

        if (nb_iter == 10000)
            break;
    }
    *vals = diag(*vals);
}
/**
 * @brief Basic QR algorithm  
 * 
 * @param M The data
 * @param vals Store the matrix containing the eigenvalues
 * @param vects Store the matrix containing the eigenvectors
 */
void eigen_decomposition(const Matrix &M, Matrix *vals, Matrix *vects)
{
    auto n = M.dim().first;
    auto m = M.dim().second;
    assert(n == m && "matrix is not square");
    *vals = M;
    *vects = identity_mat(m);
    auto nb_iter = 0;
    while (!is_triangular(*vals))
    {
        Matrix q, r;
        QR_factorization(*vals, &q, &r);
        *vals = r * q;
        *vects = *vects * q;
        nb_iter++;
        // stop when the matrix is not converging
        //most real eigenvalues matrices converge before 10000 iteration
        if (nb_iter == 10000)
            break;
    }
    *vals = diag(*vals);
}

//TODO: LU decomposition

// Singular Value decompostion
/**
 * @brief Compute the singular value decomposition of a matrix 
 * 
 * @param A The matrix to decompose
 * @param U The matrix containing the left singular vectors
 * @param Z The matrix containing the singular values
 * @param V The matrix containing the right singular values
 */

void SVD(const Matrix &A, Matrix *U, Matrix *Z, Matrix *V)
{
    // SVD: A = UZV*
    //U = eigvect(AAt), V = eigvect(AtA)
    Matrix AAt, AtA;
    AtA = A.T() * A;
    eigen_decomposition(AtA, Z, V);
    auto z = Z->dim().second;
    double a = (*Z)[0];
    bool ord = true;

    for (size_t i = 0; i != z; i++)
    {
        auto val_ii = (*Z)(0, i);
        Z->operator()(0, i) = sqrt(abs(val_ii));

        // check if Z is ordered
        if (a > val_ii)
        {
            a = val_ii;
        }
        else
            ord = false;
    }

    //check if the singular values are ordered,
    // if not, order the singular values (rare cases)
    if (ord == false)
    {
        // map indexes to their values
        map<double, size_t> ind_vals;
        for (size_t i = 0; i != z; i++)
        {
            auto value = (*Z)(0, i);
            ind_vals[value] = i;
        }

        // sort Z
        Z->sort_matrix("descending");

        // sort the indexes
        vector<size_t> indexes;
        indexes.reserve(z);
        for (size_t i = 0; i != z; i++)
        {
            indexes.emplace_back(ind_vals[(*Z)(0, i)]);
        }

        //sort V
        V->reorder_columns(indexes);
    }

    // Compute U
    *U = (A * (*V)) / (*Z);
    *V = V->T();
}

/**
 * @brief Compute the principal components of a Matrix
 * 
 * @param X The Data.
 * @param C The matrix containing the components
 * @param S The matrix containing the singular values
 * @param nb_comp The number of components
 */

void PCA(const Matrix &X, Matrix *C, Matrix *S, Matrix *Var, const size_t nb_comp = 2)
{
    Matrix U, Z, V;
    // U : left singular vectors
    // V : right singular vectors
    SVD(X, &U, &Z, &V);
    // compute explained variance
    *Var = elm2pow_n(Z, 2);
    const double s = sum(*Var);
    auto a = (*Var)[0];
    bool ord = true;
    for (auto i = 0; i != Var->nb_elm(); i++)
    {
        if (a < (*Var)[i])
            ord = true;
        (*Var)[i] = (*Var)[i] / s;
    }

    // Get the score
    *S = Z;

    // Rearrange Z to a diagonal matrix
    Z = diag_matrix(Z);

    // Make U cols and Z rows to same size, if they are not
    auto m = U.dim().second;
    auto n = Z.dim().first;
    vector_d vect_zeros;
    if (m > n)
    {
        vect_zeros.resize(m);
        for (size_t i = 0; i != m - n; i++)
            Z.add_row(vect_zeros, n);
    }
    else if (m < n)
    {
        vect_zeros.resize(n);
        for (size_t i = 0; i != m - n; i++)
            U.add_column(vect_zeros, n);
    }

    // Compute components
    *C = U * Z;

    // select desired number of components
    if (nb_comp < C->dim().second)
    {
        *C = C->sub_matrix(range(0, C->get_nb_rows() - 1), range(0, nb_comp - 1));
        *Var = Var->sub_matrix(range(0, 0), range(0, nb_comp - 1));
        *S = S->sub_matrix(range(0, 0), range(0, nb_comp - 1));
    }
}

/* OTHER FUNCTIONS */
void print(const double d)
{
    cout << d << endl;
}

void print(const string d)
{
    cout << d << endl;
}

void Matrix::print(const string s = "") const
{

    cout << s << "\n";
    for (size_t i = 1; i < nb_row * nb_col + 1; i++)
    {
        auto ij = mat_data[i - 1];
        cout.setf(ios::fixed);
        //      cout.setf(ios::showpoint);

        if (abs(ij) < 1e-10)
        {
            cout << setprecision(0);
            ij = 0;
            cout << setw(10) << ij << " ";
        }
        else if ((abs(ij) - int(abs(ij))) <= 0.000001)
        {
            cout << setprecision(0);
            cout << setw(10) << ij << " ";
        }
        else
        {
            cout << setprecision(2);
            cout << setw(10) << ij << " ";
        }

        if (i % nb_col == 0)
            cout << endl;
    }
}

void print(const Matrix &M)
{
    cout << endl;
    M.print(" ");
}

void Matrix::head(const size_t nrow = 5) const
{
    Matrix P;
    P = (*this).sub_matrix(range(0, nrow - 1), range(0, nb_col - 1));
    P.print("5 first rows:");
}

void Matrix::to_csv(const string filename) const
{
    ofstream ofile;
    ofile.open(filename);
    size_t cnt = 0;
    for (size_t i = 1; i != nb_col * nb_row + 1; i++)
    {
        ofile << mat_data[i - 1] << ",";
        long pos = ofile.tellp();
        cnt += 1;
        if (cnt % nb_col == 0)
        {
            ofile.seekp(pos - 1);
            ofile << "";
            ofile << "\n";
        }
    }
    ofile.close();
}
void save2csv(const Matrix &M, const string filename)
{
    M.to_csv(filename);
}
//TODO: Create templates;
//TODO: multithreading
