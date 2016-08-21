#include "SquareMatrix.h"

template <size_t n>
void SquareMatrix<n>::toIdentity()
{
    this->clear();
    for (size_t i = 0; i != n; ++i)
        (*this)(i, i) = 1;
}

template <size_t n>
bool SquareMatrix<n>::GetInverse(SquareMatrix<n> &invertedMatrix) const
{
    // SHOULD CHECK IF MATRIX IS SINGULAR

    /// future inverted matrix should start from identity
    invertedMatrix.toIdentity();

    /// create auxiliary matrix
    SquareMatrix<n> A = *this;

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
    } while (row-- != 0);

    return true;
}
