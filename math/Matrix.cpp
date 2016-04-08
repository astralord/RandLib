#include "Matrix.h"

template <size_t n, size_t m>
Matrix<n, m>::Matrix(double initialValue) :
    data(n * m, initialValue)
{}

template <size_t n, size_t m>
Matrix<n, m>::Matrix(const Matrix<n, m> & other) :
    data(other.data)
{}

template <size_t n, size_t m>
Matrix<n, m> & Matrix<n, m>::operator=(const Matrix<n, m> & other)
{
    if (this != & other)
    {
        data = other.data;
    }
    return *this;
}

template <size_t n, size_t m>
double &Matrix<n, m>::operator()(const size_t i, const size_t j)
{
    return data[i * m + j];
}

template <size_t n, size_t m>
double Matrix<n, m>::operator()(const size_t i, const size_t j) const
{
    return data[i * m + j];
}

template <size_t n, size_t m>
Matrix<n,m> &Matrix<n, m>::operator+=(const Matrix<n,m> &right)
{
    for (size_t i = 0; i != data.size(); ++i)
        data[i] += right.data[i];
    return *this;
}

template <size_t n, size_t m>
bool Matrix<n, m>::getSquare(Matrix<n, n> &squaredMatrix) const
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
bool Matrix<n, m>::getTransposed(Matrix<m, n> &transposedMatrix) const
{
    for (size_t i = 0; i != n; ++i)
        for (size_t j = 0; j != m; ++j)
            transposedMatrix(j, i) = (*this)(i, j);
    return true;
}

template <size_t n, size_t m>
void Matrix<n, m>::toIdentity()
{
    if (n != m) // should be part of squared
        return;
    this->clear();
    for (size_t i = 0; i != n; ++i)
        (*this)(i, i) = 1;
}

template <size_t n, size_t m>
bool Matrix<n, m>::getInverse(Matrix<m, n> &invertedMatrix) const
{
    if (n != m) // should be part of squared
        return false;

    // SHOULD CHECK IF MATRIX IS SINGULAR

    if (n == 1)
    {
        invertedMatrix(0, 0) = 1.0 / (*this)(0, 0);
        return true;
    }

    /// future inverted matrix should start from identity
    invertedMatrix.toIdentity();

    /// create auxiliary matrix
    Matrix<n, m> A = *this;

    size_t row = 0;
    while (row != n) {
        /// find first non-zero element
        size_t currentRow = row;
        while (A(currentRow, row) == 0 && currentRow < n)
            ++currentRow;
        if (currentRow == n)
            return false; /// matrix is singular

        /// swap rows if needed
        if (currentRow != row)
        {
            for (size_t i = row; i != n; ++i)
                std::swap(A(row, i), A(currentRow, i));
        }

        /// divide first row on first element
        double firstInv = 1.0 / A(row, row);
        A(row, row) = 1.0;
        invertedMatrix(row, row) = firstInv;
        for (size_t i = row + 1; i != n; ++i) {
            A(row, i) *= firstInv;
        }
        for (size_t i = 0; i != row; ++i) {
            invertedMatrix(row, i) *= firstInv;
        }

        /// subtract first row from others
        for (size_t i = row + 1; i != n; ++i) {
            double firstElement = A(i, row);

            for (size_t j = 0; j != n; ++j) {
                invertedMatrix(i, j) -= firstElement * invertedMatrix(row, j);
            }

            A(i, row) = 0.0;
            for (size_t j = row + 1; j != n; ++j) {
                A(i, j) -= firstElement * A(row, j);
            }
        }

        ++row;
    }

    /// go back
    row = n - 2;
    do {
        for (size_t i = row + 1; i != n; ++i) {
            double coef = A(row, i);
            for (size_t j = 0; j != n; ++j)
                invertedMatrix(row, j) -= invertedMatrix(i, j) * coef;
        }
    } while ((row--) != 0);

    return true;
}


template class Matrix<1, 1>;
template class Matrix<1, 2>;
template class Matrix<1, 3>;
template class Matrix<1, 4>;
template class Matrix<1, 5>;

template class Matrix<2, 1>;
template class Matrix<2, 2>;
template class Matrix<2, 3>;
template class Matrix<2, 4>;
template class Matrix<2, 5>;

template class Matrix<3, 1>;
template class Matrix<3, 2>;
template class Matrix<3, 3>;
template class Matrix<3, 4>;
template class Matrix<3, 5>;

template class Matrix<4, 1>;
template class Matrix<4, 2>;
template class Matrix<4, 3>;
template class Matrix<4, 4>;
template class Matrix<4, 5>;

template class Matrix<5, 1>;
template class Matrix<5, 2>;
template class Matrix<5, 3>;
template class Matrix<5, 4>;
template class Matrix<5, 5>;
