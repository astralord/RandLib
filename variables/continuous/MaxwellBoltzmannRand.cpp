#include "MaxwellBoltzmannRand.h"

MaxwellBoltzmannRand::MaxwellBoltzmannRand(double scale)
{
    setScale(scale);
    C.setDegree(3);
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
