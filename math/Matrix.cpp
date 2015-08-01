#include "Matrix.h"

Matrix::Matrix(unsigned height, unsigned width)
{
    n = height;
    m = width;
    data.reserve(n * m);
}

