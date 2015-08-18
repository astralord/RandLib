#include "ChiSquaredRand.h"

ChiSquaredRand::ChiSquaredRand(int degree)
{
    setDegree(degree);
}

void ChiSquaredRand::setName()
{
    nameStr = "Chi-squared(" + toStringWithPrecision(getDegree()) + ")";
}

void ChiSquaredRand::setDegree(int degree)
{
    k = std::max(degree, 1);
    if (k & 1)
    {
        pdfCoef = M_1_SQRTPI / RandMath::doubleFactorial(k - 2);
        cdfCoef = pdfCoef;
        cdfCoef *= (1 << ((k - 1) >> 1));
        pdfCoef *= M_SQRT1_2;
    }
    else
    {
        int k_2 = k >> 1;
        pdfCoef = 1.0 / RandMath::factorial(k_2 - 1);
        cdfCoef = pdfCoef;
        pdfCoef /= (1 << k_2);
    }

    GammaRand::setParameters(.5 * k, 2);
    setName();
}

double ChiSquaredRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double y = std::pow(x, .5 * k - 1);
    y *= std::exp(-.5 * x);
    return pdfCoef * y;
}

double ChiSquaredRand::F(double x) const
{
    return cdfCoef * RandMath::lowerIncGamma(.5 * k, .5 * x);
}
