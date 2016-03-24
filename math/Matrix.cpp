#include "Matrix.h"

Matrix::Matrix(const size_t height, const size_t width, const double initial_value) :
    n(height),
    m(width),
    data(n * m, initial_value)
{}

Matrix::Matrix(const Matrix & other) :
    n(other.n),
    m(other.m),
    data(other.data)
{}

Matrix & Matrix::operator=(const Matrix & other)
{
    if (this != & other)
    {
        n = other.n;
        m = other.m;
        data = other.data;
    }
    return *this;
}

double Matrix::operator()(const size_t i, const size_t j) const
{
    return data[i * n + j];
}

double & Matrix::operator()(const size_t i, const size_t j)
{
    return data[i * n + j];
}
