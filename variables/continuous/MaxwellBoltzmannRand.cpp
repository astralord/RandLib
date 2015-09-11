#include "MaxwellBoltzmannRand.h"

MaxwellBoltzmannRand::MaxwellBoltzmannRand(double scale)
{
    setScale(scale);
    C.setDegree(3);
}

std::string MaxwellBoltzmannRand::name()
{
    return "Maxwell-Boltzmann(" + toStringWithPrecision(getScale()) + ")";
}

void MaxwellBoltzmannRand::setScale(double scale)
{
    a = std::max(scale, MIN_POSITIVE);
}

double MaxwellBoltzmannRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double x2 = x * x;
    double a2 = a * a;
    double y = std::exp(-.5 * x2 / a2);
    return M_SQRT2 * M_1_SQRTPI * x2 * y / (a2 * a);
}

double MaxwellBoltzmannRand::F(double x) const
{
    if (x <= 0)
        return 0;
    double xAdj = M_SQRT1_2 * x / a;
    double y = std::exp(-xAdj * xAdj);
    y *= M_SQRT2 * M_1_SQRTPI * x / a;
    return std::erf(xAdj) - y;
}

double MaxwellBoltzmannRand::variate() const
{
    return a * std::sqrt(C.variate());
}

double MaxwellBoltzmannRand::E() const
{
    return 2 * M_1_SQRTPI * M_SQRT2 * a;
}

double MaxwellBoltzmannRand::Var() const
{
    return a * a * (3 - 8.0 * M_1_PI);
}

double MaxwellBoltzmannRand::Mode() const
{
    return M_SQRT2 * a;
}

double MaxwellBoltzmannRand::Skewness() const
{
    double skewness = 3 * M_PI - 8;
    skewness = 2.0 / skewness;
    skewness *= std::sqrt(skewness);
    return (16 - 5 * M_PI) * skewness;
}

double MaxwellBoltzmannRand::ExcessKurtosis() const
{
    double numerator = 40 - 3 * M_PI;
    numerator *= M_PI;
    numerator -= 96;
    double denominator = 3 * M_PI - 8;
    denominator *= denominator;
    return 4 * numerator / denominator;
}
