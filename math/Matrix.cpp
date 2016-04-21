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
        data = other.data;
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
double &Matrix<n, m>::operator()(const size_t i)
{
    return data[i];
}

template <size_t n, size_t m>
double Matrix<n, m>::operator()(const size_t i) const
{
    return data[i];
}

template <size_t n, size_t m>
Matrix<n,m> &Matrix<n, m>::operator+=(const Matrix<n,m> &right)
{
    for (size_t i = 0; i != data.size(); ++i)
        data[i] += right.data[i];
    return *this;
}

template <size_t n, size_t m>
void Matrix<n, m>::getSquare(Matrix<n, n> &squaredMatrix) const
{
    for (size_t i = 0; i != n; ++i) {
        for (size_t j = i; j != n; ++j) {
            squaredMatrix(i, j) = 0.0;
            for (size_t k = 0; k != m; ++k)
                squaredMatrix(i, j) += (*this)(i, k) * (*this)(k, j);
            squaredMatrix(j, i) = squaredMatrix(i, j);
        }
    }
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

/// define Matrix(N, M + [1, 10])
#define DEFINE_VECTOR_DECADE(N, M) \
    template class Matrix<N, M + 1>; \
    template class Matrix<N, M + 2>; \
    template class Matrix<N, M + 3>; \
    template class Matrix<N, M + 4>; \
    template class Matrix<N, M + 5>; \
    template class Matrix<N, M + 6>; \
    template class Matrix<N, M + 7>; \
    template class Matrix<N, M + 8>; \
    template class Matrix<N, M + 9>; \
    template class Matrix<N, M + 10>;

/// define Matrix(N + [1, 10], M + [1, 10])
#define DEFINE_MATRIX_DECADE(N, M) \
    DEFINE_VECTOR_DECADE(N + 1, M) \
    DEFINE_VECTOR_DECADE(N + 2, M) \
    DEFINE_VECTOR_DECADE(N + 3, M) \
    DEFINE_VECTOR_DECADE(N + 4, M) \
    DEFINE_VECTOR_DECADE(N + 5, M) \
    DEFINE_VECTOR_DECADE(N + 6, M) \
    DEFINE_VECTOR_DECADE(N + 7, M) \
    DEFINE_VECTOR_DECADE(N + 8, M) \
    DEFINE_VECTOR_DECADE(N + 9, M) \
    DEFINE_VECTOR_DECADE(N + 10, M)

DEFINE_MATRIX_DECADE(0, 0)
DEFINE_MATRIX_DECADE(0, 10)
DEFINE_MATRIX_DECADE(10, 0)
DEFINE_MATRIX_DECADE(10, 10)


