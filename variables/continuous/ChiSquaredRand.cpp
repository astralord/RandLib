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
    GammaRand::setParameters(.5 * std::max(degree, 1), 2);
}
