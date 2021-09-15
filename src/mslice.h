#ifndef MSLICE_H
#define MSLICE_H
//#include <iostream>
#include "../include/matrix.h"
//struct range;

//class Matrix; 
struct Mslice
{
    Matrix *ptr;
    range row;
    range col;
    Mslice() = default;
    Mslice(const range r, const range c, Matrix *p);
    virtual ~Mslice();
    Mslice(const Mslice &other) = default;
    Matrix &operator=(const initializer_list<vector<double>> li);
    Matrix &operator=(const Matrix &);
};

#endif // MSLICE_H
