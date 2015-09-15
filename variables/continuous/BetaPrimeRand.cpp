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
    double rv = std::pow(x, alpha - 1);
    rv *= std::pow(1 + x, -alpha - beta);
    return pdfCoef * rv;
}

double BetaPrimeRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return BetaRand::F(x / (1.0 + x));
}

double BetaPrimeRand::variate() const
{
    double x = BetaRand::variate();
    return x / (1.0 - x);
}

void BetaPrimeRand::sample(QVector<double> &outputData) const
{
    BetaRand::sample(outputData);
    for (double &var : outputData)
        var = var / (1.0 - var);
}

double BetaPrimeRand::Mean() const
{
    return (beta > 1) ? alpha / (beta - 1) : INFINITY;
}

double BetaPrimeRand::Variance() const
{
    if (beta <= 2)
        return INFINITY;
    double betaAdj = beta - 1;
    double numerator = alpha * (alpha + betaAdj);
    double denominator = (betaAdj - 1) * betaAdj * betaAdj;
    return numerator / denominator;
}

double BetaPrimeRand::Median() const
{
    double betaMedian = BetaRand::Median();
    qDebug() << betaMedian;
    return betaMedian / (1.0 - betaMedian);
}

double BetaPrimeRand::Mode() const
{
    return (alpha < 1) ? 0 : (alpha - 1) / (beta + 1);
}

double BetaPrimeRand::Skewness() const
{
    if (beta <= 3)
        return INFINITY;
    double aux = alpha + beta - 1;
    double skewness = (beta - 2) / (alpha * aux);
    skewness = std::sqrt(skewness);
    aux += alpha;
    aux += aux;
    return aux * skewness / (beta - 3);
}

double BetaPrimeRand::ExcessKurtosis() const
{
    if (beta <= 4)
        return INFINITY;
    double betam1 = beta - 1;
    double numerator = betam1 * betam1 * (beta - 2) / (alpha * (alpha + betam1));
    numerator += 5 * beta - 11;
    double denominator = (beta - 3) * (beta - 4);
    return 6 * numerator / denominator;
}
