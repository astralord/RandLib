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

    Y.setParameters(m, w / m);

    sigma = m / w;
    logGammaM = Y.getLogGammaFunction();
    pdfCoef = M_LN2 + m * std::log(m / w) - logGammaM;
}

double NakagamiRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double y = (m + m - 1) * std::log(x);
    y -= sigma * x * x;
    return std::exp(pdfCoef + y);
}

double NakagamiRand::F(double x) const
{
    if (x <= 0)
        return 0;
    double y = RandMath::logLowerIncGamma(m, sigma * x * x);
    return std::exp(y - logGammaM);
}

double NakagamiRand::variate() const
{
    return std::sqrt(Y.variate());
}

double NakagamiRand::Mean() const
{
    double y = std::lgamma(m + 0.5);
    y -= logGammaM;
    y -= 0.5 * (pdfCoef - M_LN2 + logGammaM) / m;
    return std::exp(y);
}

double NakagamiRand::Variance() const
{
    double y = std::lgamma(m + 0.5);
    y -= logGammaM;
    y = std::exp(y + y);
    return w * (1 - y / m);
}

double NakagamiRand::Mode() const
{
    double mode = 0.5 * w / m;
    return std::sqrt(w - mode);
}
