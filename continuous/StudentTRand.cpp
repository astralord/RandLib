#include "StudentTRand.h"

StudentTRand::StudentTRand(int degree) :
    X(0, 1),
    Y(degree)
{
    setDegree(degree);
}

void StudentTRand::setDegree(double degree)
{
    v = std::max(degree, MIN_POSITIVE);
    Y.setDegree(v);

    pdfCoef = std::tgamma(.5 * (v + 1));
    pdfCoef /= (qSqrt(v * M_PI) * std::tgamma(.5 * v));
}

double StudentTRand::pdf(double x) const
{
    double y = 1 + x * x / v;
    y = std::pow(y, -.5 * (v + 1));
    return pdfCoef * y;
}

double StudentTRand::cdf(double x) const
{
    return x;
}

double StudentTRand::value()
{
    //v = 1 -> cauchy
    //v = inf -> normal
    return X.value() * qSqrt((double)v / Y.value());
}
