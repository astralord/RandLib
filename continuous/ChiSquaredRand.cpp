#include "ChiSquaredRand.h"

ChiSquaredRand::ChiSquaredRand(int k) :
    X(0, 1)
{
    setDegree(k);
}

void ChiSquaredRand::setDegree(int degree)
{
    k = std::max(degree, 1);
    if (k % 2 == 0)
    {
        size_t k_2 = .5 * k;
        pdfCoef = 1.0 / RandMath::fastFactorial(k_2 - 1);
        cdfCoef = pdfCoef;
        pdfCoef /= (1 << k_2);
    }
    else
    {
        pdfCoef = 1.0 / RandMath::doubleFactorial(k - 2);
        cdfCoef = pdfCoef;
        cdfCoef *= M_1_SQRTPI;
        cdfCoef *= (1 << (size_t)(.5 * (k - 1)));
        pdfCoef *= M_1_SQRT2PI;
    }
}

double ChiSquaredRand::pdf(double x) const
{
    if (x <= 0)
        return 0;
    double y = std::pow(x, .5 * k - 1);
    y *= std::exp(-.5 * x);
    return pdfCoef * y;
}

double ChiSquaredRand::cdf(double x) const
{
    return cdfCoef * RandMath::lowerIncGamma(.5 * k, .5 * x);
}

double ChiSquaredRand::value()
{
    // TODO: Check if this can be replaced with gamma distribution generator
    // and maybe it should if k is too big
    // tip: Y ~ Gamma(k / 2, 2)

    double rv = 0;
    for (int i = 0; i != k; ++i)
    {
        double x = X.value();
        x *= x;
        rv += x;
    }
    return rv;
}
