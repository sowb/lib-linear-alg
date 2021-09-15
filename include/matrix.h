#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm> //generate
#include <cmath>
#include <map>
#include <iomanip> // std::setprecision
#include <cassert> /* assert */
#include <random>
#include <memory>
#include <functional>

using namespace std;

typedef vector<double> vector_d;
typedef _Bind<double (*(double, double))(double, double)> bind_d;
struct range
{
    size_t start;
    size_t end;
    size_t len;
    range(size_t x, size_t y) : start(x), end(y + 1), len(end - start) {}
};

struct Mslice;

class Matrix
{
    vector_d mat_data;
    size_t nb_row = 0;
    size_t nb_col = 0;
    inline void split_str(string &str, string delimiter);

public:
    //Constructors
    Matrix() = default;
    Matrix(Matrix &&) = default;
    Matrix &operator=(Matrix &&) = default;
    Matrix(const string, const string, const bool);       // Load CVS files
    Matrix(const size_t n, const size_t m, bind_d); //FIXME: can be optimised
    Matrix(const size_t n, const size_t m, double fill);
    Matrix(initializer_list<std::vector<double>> li);
    Matrix(const Matrix &);
    virtual ~Matrix() {}

    // Operators Overloading
    double &operator()(const size_t n, const size_t m);
    double operator()(size_t n, size_t m) const;
    Mslice operator()(const range row, const range col);
    double &operator[](const size_t i);
    double operator[](const size_t) const;
    Matrix &operator=(const Matrix &);
    friend bool operator==(const Matrix &lmat, const Matrix &rmat);
    friend bool operator!=(const Matrix &lmat, const Matrix &rmat);
    Matrix &operator+=(const Matrix &);
    Matrix &operator-=(const Matrix &);
    Matrix &operator*=(const double val);

    // Members for special matrices and conversions
    Matrix Id(const size_t n);
    //void reserve_m(size_t n) { mat_data.reserve(n); }
    //friend Matrix transpose(const Matrix &);
    //Matrix T();
    Matrix T() const;
    friend Matrix as_matrix(vector<double>, size_t);
    inline friend vector_d as_vector(const Matrix &);
    friend Matrix elm2pow_n(const Matrix &, double);

    // Members and friends functions to slice, subscript a matrix
    inline size_t nb_elm() const { return nb_row * nb_col; }
    inline size_t idx(size_t i, size_t j, size_t stride) { return i * stride + j; } //TODO: add inline
    inline size_t idx(size_t i, size_t j, size_t stride) const { return i * stride + j; }
    Matrix sub_matrix(const range row, const range col) const;
    vector_d column(size_t m, size_t starting_row);
    vector_d column(size_t m, size_t starting_row) const;
    vector_d row(size_t n);
    vector_d row(size_t n) const;
    void add_row(const vector<double>, const size_t pos);
    void add_column(const vector<double>, const size_t pos);
    void remove_column(const size_t);
    void shape();
    pair<size_t, size_t> dim() { return make_pair(nb_row, nb_col); }
    pair<size_t, size_t> dim() const { return make_pair(nb_row, nb_col); }
    void reshape(const size_t n, const size_t m);
    const size_t get_nb_rows() const { return nb_row; }
    const size_t get_nb_columns() const { return nb_col; }
    void sort_matrix(string);
    void reorder_columns(vector<size_t>);
    void print(const string s) const;
    void head(const size_t) const;
    void to_csv(const string) const;
    //Member for computing statistics
    inline double sum() const { return accumulate(mat_data.begin(), mat_data.end(), 0.0); }
    inline double avg() const { return this->sum() / mat_data.size(); }
};

/* Nom Member Declaration; */
void print(const double);
void print(const string);
void print(const Matrix &);
void save2csv(const Matrix &, const string);

//NON MEMBER OPERATORS
bool operator==(const Matrix &, const Matrix &);
bool operator!=(const Matrix &, const Matrix &);
Matrix operator+(const Matrix &, const Matrix &);
Matrix operator+(const Matrix &, const double);
Matrix operator+(const double, const Matrix &);

Matrix operator-(const Matrix &, const Matrix &);
Matrix operator-(const Matrix &, const double);
Matrix operator-(const double, const Matrix &);

Matrix operator*(const Matrix &, const Matrix &);
Matrix operator*(const Matrix &, const double);
Matrix operator*(const double, const Matrix &);

Matrix operator/(const Matrix &, const double);
Matrix operator/(const Matrix &, const Matrix &);

//MATRIX OPERATIONS
bool desc(double, double);
Matrix identity_mat(const size_t);
Matrix diag_matrix(const Matrix &);
Matrix diag(const Matrix &);
Matrix transpose(const Matrix &M);
double sgn(double val);
double norm(const vector<double> &);
double norm(const Matrix &);
bool is_triangular(const Matrix &T);

//STAT
Matrix elm2pow_n(const Matrix &M, double);
double sum(const Matrix &);
Matrix tile(Matrix, const size_t);
double avg(vector<double>);
Matrix mean(const Matrix &, const string);
Matrix stdev(const Matrix &, const string);
Matrix scale(const Matrix &, const string);
double uniform_dist(double, double);
double normal_dist(double, double);
bind_d unif(double, double);
bind_d rnorm(double, double);
//ALGORITHMS
Matrix householder(Matrix);
Matrix hessenberg(const Matrix &);
void QR_factorization(Matrix, Matrix *, Matrix *);
void eigen_decomposition(const Matrix &, Matrix *, Matrix *);
double wilkinson_shift(const double, const double, const double c);
void QR_algorithm_w_shift(const Matrix &, Matrix *, Matrix *);
void SVD(const Matrix &, Matrix *, Matrix *, Matrix *);
void PCA(const Matrix &, Matrix *, Matrix *, Matrix *, const size_t);

#endif