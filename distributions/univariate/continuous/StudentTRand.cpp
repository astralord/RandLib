#include "StudentTRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"

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

    vp1Half = 0.5 * (v + 1);
    pdfCoef = std::lgamma(vp1Half);
    pdfCoef -= 0.5 * std::log(M_PI * v);
    pdfCoef -= Y.getLogGammaFunction();
}

double StudentTRand::f(double x) const
{
    double y = 1 + x * x / v;
    y = -vp1Half * std::log(y);
    return std::exp(pdfCoef + y);
}

double StudentTRand::F(double x) const
{
    if (v == 1)
        return 0.5 + std::atan(x) * M_1_PI;
    if (v == 2)
        return 0.5 * (1.0 + x / std::sqrt(2 + x * x));
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

void StudentTRand::sample(std::vector<double> &outputData) const
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

double StudentTRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
         return NAN;
    if (p == 0)
        return -INFINITY;
    if (p == 1)
        return INFINITY;

    double temp = p - 0.5;
    if (v == 1)
        return std::tan(M_PI * temp);
    double alpha = 4 * p * (1.0 - p);
    if (v == 2)
        return 2.0 * temp * std::sqrt(2.0 / alpha);
    if (v == 4)
    {
        double sqrtAlpha = std::sqrt(alpha);
        double q = std::cos(std::acos(sqrtAlpha) / 3.0) / sqrtAlpha;
        return 2 * RandMath::sign(temp) * std::sqrt(q - 1);
    }
    return ContinuousDistribution::Quantile(p);
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
