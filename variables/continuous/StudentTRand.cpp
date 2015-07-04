#include "StudentTRand.h"

StudentTRand::StudentTRand(int degree)
{
    setDegree(degree);
}

void StudentTRand::setDegree(double degree)
{
    v = std::max(degree, MIN_POSITIVE);
    Y.setDegree(static_cast<int>(v));

    pdfCoef = std::tgamma(.5 * (v + 1));
    pdfCoef /= (std::sqrt(v * M_PI) * std::tgamma(.5 * v));
}

double StudentTRand::f(double x) const
{
    double y = 1 + x * x / v;
    y = std::pow(y, -.5 * (v + 1));
    return pdfCoef * y;
}

double StudentTRand::F(double x) const
{
    return x;
}

double StudentTRand::value()
{
    //v = 1 -> cauchy
    //v = inf -> normal
    return X.value() * std::sqrt(v / Y.value());
}