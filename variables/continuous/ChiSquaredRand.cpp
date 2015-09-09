#include "ChiSquaredRand.h"

ChiSquaredRand::ChiSquaredRand(int degree)
{
    setDegree(degree);
}

std::string ChiSquaredRand::name()
{
    return "Chi-squared(" + toStringWithPrecision(getDegree()) + ")";
}

void ChiSquaredRand::setDegree(size_t degree)
{
    double shape = (degree < 1) ? 0.5 : 0.5 * degree;
    GammaRand::setParameters(shape, 2);
}
