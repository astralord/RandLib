#include "ChiSquaredRand.h"

ChiSquaredRand::ChiSquaredRand(int degree)
{
    setDegree(degree);
}

std::string ChiSquaredRand::name()
{
    return "Chi-squared(" + toStringWithPrecision(getDegree()) + ")";
}

void ChiSquaredRand::setDegree(int degree)
{
    k = std::max(degree, 1);
    GammaRand::setParameters(.5 * k, 2);
}
