#include "NakagamiRand.h"

NakagamiRand::NakagamiRand(double shape, double spread)
{
    setParameters(shape, spread);
}

std::string NakagamiRand::name()
{
    return "Nakagami(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getSpread()) + ")";
}

void NakagamiRand::setParameters(double shape, double spread)
{
    m = std::max(shape, 0.5);
    w = spread;
    if (w <= 0)
        w = 1.0;

    sigma = m / w;
    gammaMInv = 1.0 / std::tgamma(m);
    pdfCoef = 2.0 * std::pow(m / w, m) * gammaMInv;

    Y.setParameters(m, w / m);
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

double NakagamiRand::Mean() const
{
    return std::tgamma(m + 0.5) * gammaMInv * std::sqrt(w / m);
}

double NakagamiRand::Variance() const
{
    double res = std::tgamma(m + 0.5) * gammaMInv;
    res *= res;
    return w * (1 - res / m);
}

double NakagamiRand::Mode() const
{
    double mode = 0.5 * w / m;
    return std::sqrt(w - mode);
}
