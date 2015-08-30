#include "BetaPrimeRand.h"

BetaPrimeRand::BetaPrimeRand(double shape1, double shape2)
    : BetaRand(shape1, shape2)
{
}

std::string BetaPrimeRand::name()
{
    return "Beta Prime(" + toStringWithPrecision(getAlpha()) + ", " + toStringWithPrecision(getBeta()) + ")";
}

double BetaPrimeRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double alpha = X.getShape();
    double beta = Y.getShape();
    double rv = std::pow(x, alpha - 1);
    rv *= std::pow(1 + x, -alpha - beta);
    return pdfCoef * rv;
}

double BetaPrimeRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return RandMath::integral([this] (double t)
    {
        return f(t);
    },
    0, x);
}

double BetaPrimeRand::variate() const
{
    double x = BetaRand::variate();
    return x / (1.0 - x);
}

void BetaPrimeRand::sample(QVector<double> &outputData)
{
    BetaRand::sample(outputData);
    for (double &var : outputData)
        var = var / (1.0 - var);
}
