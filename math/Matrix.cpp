#include "Matrix.h"

Matrix::Matrix(const size_t height, const size_t width, const double initial_value) :
    n(height),
    m(width),
    data(n * m, initial_value),
    isTransposed(false)
{}

Matrix::Matrix(const Matrix & other) :
    n(other.n),
    m(other.m),
    data(other.data),
    isTransposed(false)
{}

Matrix & Matrix::operator=(const Matrix & other)
{
    if (this != & other)
    {
        n = other.n;
        m = other.m;
        data = other.data;
        isTransposed = other.isTransposed;
    }
    return *this;
}

double Matrix::operator()(const size_t i, const size_t j) const
{
    return (isTransposed) ? data[j * n + i] : data[i * m + j];
}

double & Matrix::operator()(const size_t i, const size_t j)
{
    return (isTransposed) ? data[j * n + i] : data[i * m + j];
}

Matrix &Matrix::operator+=(const Matrix &right)
{
    if (n != right.height() || m != right.width())
        return *this; // we should throw exception here
    for (size_t i = 0; i != data.size(); ++i)
        data[i] += right.data[i];
    return *this;
}

bool Matrix::getSquare(SquareMatrix &squaredMatrix) const
{
    if (squaredMatrix.size() != n)
        return false;
    for (size_t i = 0; i != n; ++i) {
        for (size_t j = i; j != n; ++j) {
            squaredMatrix(i, j) = 0.0;
            for (size_t k = 0; k != m; ++k)
                squaredMatrix(i, j) += (*this)(i, k) * (*this)(k, j);
            squaredMatrix(j, i) = squaredMatrix(i, j);
        }
    }
    return true;
}

void Matrix::clear()
{
    std::fill(data.begin(), data.end(), 0);
}

void Matrix::transpose()
{
    std::swap(n, m);
    isTransposed = !isTransposed;
}


/// VECTOR
Vector::Vector(const size_t height, const double initial_value) :
    Matrix(height, 1, initial_value)
{
}

Vector &Vector::operator+=(const Vector &right)
{
    if (n != right.height())
        return *this; // we should throw exception here
    for (size_t i = 0; i != data.size(); ++i)
        data[i] += right.data[i];
    return *this;
}

/// SQUARE MATRIX
SquareMatrix::SquareMatrix(const size_t height, const double initial_value) :
    Matrix(height, height, initial_value)
{
}

void SquareMatrix::toIdentity()
{
    clear();
    for (size_t i = 0; i != n; ++i)
        (*this)(i, i) = 1;
}

bool SquareMatrix::getInverse(SquareMatrix & invertedMatrix)
{
    if (invertedMatrix.size() != n)
        return false;

    // SHOULD CHECK IF MATRIX IS SINGULAR

    if (n == 1)
    {
        invertedMatrix(0, 0) = 1.0 / (*this)(0, 0);
        return true;
    }

    /// future inverted matrix should start from identity
    invertedMatrix.toIdentity();

    size_t row = 0;
    while (row != n) {
        /// find first non-zero element
        size_t currentRow = row;
        while ((*this)(currentRow, row) == 0 && currentRow < n)
            ++currentRow;
        if (currentRow == n)
            return false; /// matrix is singular

        /// swap rows if needed
        if (currentRow != row)
        {
            for (size_t i = row; i != n; ++i)
                std::swap((*this)(row, i), (*this)(currentRow, i));
        }

        /// divide first row on first element
        double firstInv = 1.0 / (*this)(row, row);
        (*this)(row, row) = 1.0;
        invertedMatrix(row, row) = firstInv;
        for (size_t i = row + 1; i != n; ++i) {
            (*this)(row, i) *= firstInv;
        }
        for (size_t i = 0; i != row; ++i) {
            invertedMatrix(row, i) *= firstInv;
        }

        /// subtract first row from others
        for (size_t i = row + 1; i != n; ++i) {
            double firstElement = (*this)(i, row);

            for (size_t j = 0; j != n; ++j) {
                invertedMatrix(i, j) -= firstElement * invertedMatrix(row, j);
            }

            (*this)(i, row) = 0.0;
            for (size_t j = row + 1; j != n; ++j) {
                (*this)(i, j) -= firstElement * (*this)(row, j);
            }
        }

        ++row;
    }

    /// go back
    row = n - 2;
    do {
        for (size_t i = row + 1; i != n; ++i) {
            double coef = (*this)(row, i);
            for (size_t j = 0; j != n; ++j)
                invertedMatrix(row, j) -= invertedMatrix(i, j) * coef;
        }
    } while ((row--) != 0);

    return true;
}

