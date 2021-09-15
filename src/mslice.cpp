#include "mslice.h"

Mslice::Mslice(const range r, const range c, Matrix *p) : row(r), col(c), ptr(p) {}

Mslice::~Mslice()
{
}

Matrix &Mslice::operator=(const initializer_list<vector<double>> li)
{
    vector<vector<double>> D(li);
    assert(make_pair(D.size(), D[0].size()) == make_pair(row.len, col.len));
    size_t r = 0; // row index for rhs

    for (size_t i = row.start; i != row.end; i++)
    {
        size_t c = 0;
        for (size_t j = col.start; j != col.end; j++)
        {
            (*ptr)(i, j) = D[r][c];
            c++;
        }
        r++;
    }
    return *ptr;
}

Matrix &Mslice::operator=(const Matrix &rhs)
{
    if (rhs.dim() != make_pair(row.len, col.len))
        throw runtime_error("Not same dim");

    size_t r = 0; // row index for rhs
    for (size_t i = row.start; i != row.end; i++)
    {
        size_t c = 0; // col index for rhs
        for (size_t j = col.start; j != col.end; j++)
        {
            (*ptr)(i, j) = rhs(r, c);
            c++;
        }
        r++;
    }
    return *ptr;
}
