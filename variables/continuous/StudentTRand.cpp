#include "StudentTRand.h"

StudentTRand::StudentTRand(int degree)
{
    setDegree(degree);
}

std::string StudentTRand::name()
{
    return "Student's t(" + toStringWithPrecision(getDegree()) + ")";
}

void StudentTRand::setDegree(int degree)
{
    v = std::max(degree, 1);
    Y.setDegree(v);

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
    double absY = RandMath::integral([this] (double t)
    {
        return f(t);
    },
    0, std::fabs(x));
    return (x > 0) ? 0.5 + absY : 0.5 - absY;
}

double StudentTRand::variate() const
{
    //v = inf -> normal
    if (v == 1)
        return CauchyRand::standardVariate();
    return NormalRand::standardVariate() / std::sqrt(Y.variate() / v);
}

void StudentTRand::sample(QVector<double> &outputData)
{
    //v = inf -> normal
    if (v == 1) {
        for (double &var : outputData)
            var = CauchyRand::standardVariate();
    }
    else {
        for (double &var : outputData)
            var = NormalRand::standardVariate() / std::sqrt(Y.variate() / v);
    }
}

double StudentTRand::Mean() const
{
    return (v > 1) ? 0 : NAN;
}

double StudentTRand::Variance() const
{
    if (v > 2)
        return v / (v - 2);
    return (v > 1) ? INFINITY : NAN;
}

double StudentTRand::Median() const
{
    return 0.0;
}

double StudentTRand::Mode() const
{
    return 0.0;
}

double StudentTRand::Skewness() const
{
    return (v > 3) ? 0.0 : NAN;
}

double StudentTRand::ExcessKurtosis() const
{
    if (v > 4)
        return 6 / (v - 4);
    return (v > 2) ? INFINITY : NAN;
}
