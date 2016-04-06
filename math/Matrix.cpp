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
    return data[i * m + j];
}

double & Matrix::operator()(const size_t i, const size_t j)
{
    return data[i * m + j];
}

Matrix &Matrix::operator+=(const Matrix &right)
{
    if (height() != right.height() || width() != right.width())
        return *this; // we should throw exception here
    for (size_t i = 0; i != data.size(); ++i)
        data[i] += right.data[i];
    return *this;
}
