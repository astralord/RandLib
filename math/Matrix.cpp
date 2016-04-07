#include "Matrix.h"

template <size_t n, size_t m>
Matrix<n, m>::Matrix(double initialValue) :
    data(n * m, initialValue),
    isTransposed(false)
{}

template <size_t n, size_t m>
Matrix<n, m>::Matrix(const Matrix<n, m> & other) :
    data(other.data),
    isTransposed(false)
{}

template <size_t n, size_t m>
Matrix<n, m> & Matrix<n, m>::operator=(const Matrix<n, m> & other)
{
    if (this != & other)
    {
        data = other.data;
        isTransposed = other.isTransposed;
    }
    return *this;
}

template <size_t n, size_t m>
Matrix<n,m> &Matrix<n, m>::operator+=(const Matrix<n,m> &right)
{
    for (size_t i = 0; i != data.size(); ++i)
        data[i] += right.data[i];
    return *this;
}

template <size_t n, size_t m>
bool Matrix<n, m>::getSquare(SquareMatrix<n> &squaredMatrix) const
{
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

template <size_t n, size_t m>
void Matrix<n, m>::fill(double value)
{
    std::fill(data.begin(), data.end(), value);
}

template <size_t n, size_t m>
void Matrix<n, m>::clear()
{
    fill(0);
}

template <size_t n, size_t m>
void Matrix<n, m>::transpose()
{
    size_t c = n;
    n = m;
    m = c;
    isTransposed = !isTransposed;
}

/// SQUARE MATRIX
template <size_t n>
SquareMatrix<n>::SquareMatrix(const double initial_value) :
    Matrix<n, n>(initial_value)
{
}

template <size_t n>
void SquareMatrix<n>::toIdentity()
{
    this->clear();
    for (size_t i = 0; i != n; ++i)
        (*this)(i, i) = 1;
}

template <size_t n>
bool SquareMatrix<n>::getInverse(SquareMatrix<n> &invertedMatrix) const
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

