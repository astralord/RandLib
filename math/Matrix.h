#ifndef MATRIX_H
#define MATRIX_H

#include "randlib_global.h"
#include <QVector>
#include <QDebug>

/**
 * @brief The Matrix class
 *
 */
class RANDLIBSHARED_EXPORT Matrix
{
    unsigned n, m;
    QVector<double> data;
public:
    Matrix(unsigned height, unsigned width);
    inline unsigned getHeight() { return n; }
    inline unsigned getWidth() { return m; }

    double &operator()(unsigned i, unsigned j) {
        return data[i * n + j];
    }

    // TODO:
    // addition and multiplication
};


#endif // MATRIX_H
