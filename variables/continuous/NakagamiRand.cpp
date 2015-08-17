#include "NakagamiRand.h"

NakagamiRand::NakagamiRand(double shape, double spread)
{
    setParameters(shape, spread);
}

void NakagamiRand::setName()
{
    nameStr = "Nakagami(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getSpread()) + ")";
}

void NakagamiRand::setParameters(double shape, double spread)
{
    m = std::max(shape, 0.5);
    w = std::max(spread, MIN_POSITIVE);

    sigma = m / w;
    gammaMInv = 1.0 / std::tgamma(m);
    pdfCoef = 2.0 * std::pow(m / w, m) * gammaMInv;

    Y.setParameters(m, w / m);
    setName();
}

double NakagamiRand::f(double x) const
{
    if (x <= 0)
        return 0;
    return pdfCoef * std::pow(x, 2.0 * m - 1) * std::exp(-sigma * x * x);
}

double NakagamiRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return RandMath::lowerIncGamma(m, sigma * x * x) * gammaMInv;
}

double NakagamiRand::variate() const
{
    return std::sqrt(Y.variate());
}

