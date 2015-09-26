#include "BinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/NormalRand.h"
#include "../continuous/ExponentialRand.h"
#include "GeometricRand.h"

BinomialRand::BinomialRand(int number, double probability)
{
    setParameters(number, probability);
}

std::string BinomialRand::name()
{
    return "Binomial(" + toStringWithPrecision(getNumber()) + ", " + toStringWithPrecision(getProbability()) + ")";
}

void BinomialRand::setGeneratorConstants()
{
    double minpq = std::min(p, q);
    npFloor = std::floor(n * minpq);
    if (npFloor <= generatorEdge)
        return;
    nqFloor = n - npFloor;
    pFloor = npFloor / n;
    pRes = minpq - pFloor;
    double qFloor = 1.0 - pFloor;

    /// set deltas
    double npq = npFloor * qFloor;
    double coef = 128.0 * n / M_PI;
    delta1 = coef * pFloor / (81.0 * qFloor);
    delta1 = npq * std::log(delta1);
    if (delta1 > 1.0)
        delta1 = std::sqrt(delta1);
    else
        delta1 = 1.0;
    delta2 = coef * qFloor / pFloor;
    delta2 = npq * std::log(delta2);
    if (delta2 > 1.0)
        delta2 = std::sqrt(delta2);
    else
        delta2 = 1.0;

    /// set sigmas and c
    double npqSqrt = std::sqrt(npq);
    sigma1 = npqSqrt * (1.0 + 0.25 * delta1 / npFloor);
    sigma2 = npqSqrt * (1.0 + 0.25 * delta2 / nqFloor);
    c = 2.0 * delta1 / npFloor;

    /// set a's
    a1 = 0.5 * std::exp(c) * sigma1 * M_SQRT2PI;

    a2 = 0.5 * sigma2 * M_SQRT2PI;
    a2 += a1;

    coefa3 = 0.5 * delta1 / (sigma1 * sigma1);
    a3 = 1.0 / nqFloor - coefa3;
    a3 *= delta1;
    a3 = std::exp(a3);
    a3 /= coefa3;
    a3 += a2;

    coefa4 = 0.5 * delta2 / (sigma2 * sigma2);
    a4 = std::exp(-delta2 * coefa4) / coefa4;
    a4 += a3;

    PnpInv = 1.0 / PFloor(npFloor);
}

void BinomialRand::setParameters(int number, double probability)
{
    n = std::max(number, 1);
    p = std::min(probability, 1.0);
    p = std::max(p, 0.0);
    q = 1.0 - p;
    np = n * p;
    setGeneratorConstants();
}

double BinomialRand::PFloor(int k) const
{
    if (k == n - k)
        return RandMath::binomialCoef(n, k) * std::pow(pFloor * (1.0 - pFloor), k);
    return RandMath::binomialCoef(n, k) * std::pow(pFloor, k) * std::pow(1.0 - pFloor, n - k); // TODO: storage two logs and don't calculate pow
}

double BinomialRand::P(int k) const
{
    if (k < 0 || k > n)
        return 0;
    if (k == n - k)
        return RandMath::binomialCoef(n, k) * std::pow(p * q, k);
    return RandMath::binomialCoef(n, k) * std::pow(p, k) * std::pow(q, n - k);
}

double BinomialRand::F(double x) const
{
    if (x < 0)
        return 0;
    int k = std::floor(x);
    return RandMath::regularizedBetaFun(q, n - k, 1 + k);
}

double BinomialRand::variateRejection() const
{
    bool reject = true;
    int X = 0, iter = 0;
    double Y, V;
    do {
        double U = UniformRand::variate(0, a4);
        if (U <= a1)
        {
            double N = NormalRand::standardVariate();
            Y = sigma1 * std::fabs(N);
            reject = (Y >= delta1);
            if (!reject)
            {
                double W = ExponentialRand::standardVariate();
                X = static_cast<int>(std::floor(Y));
                V = -W - 0.5 * N * N + c;
            }
        }
        else if (U <= a2)
        {
            double N = NormalRand::standardVariate();
            Y = sigma2 * std::fabs(N);
            reject = (Y >= delta2);
            if (!reject)
            {
                double W = ExponentialRand::standardVariate();
                X = static_cast<int>(std::floor(-Y));
                V = -W - 0.5 * N * N;
            }
        }
        else if (U <= a3)
        {
            double W1 = ExponentialRand::standardVariate();
            double W2 = ExponentialRand::standardVariate();
            Y = delta1 + W1 / coefa3;
            X = static_cast<int>(std::floor(Y));
            V = -W2 - coefa3 * Y + delta1 / nqFloor;
            reject = false;
        }
        else
        {
            double W1 = ExponentialRand::standardVariate();
            double W2 = ExponentialRand::standardVariate();
            Y = delta2 + W1 / coefa4;
            X = static_cast<int>(std::floor(-Y));
            V = -W2 - coefa4 * Y;
            reject = false;
        }

        X += npFloor;
        if (!reject && X >= 0 && X <= n && V <= std::log(PFloor(X) * PnpInv))
            return X;

    } while (++iter < 1e9);
    return NAN;
}

double BinomialRand::variateWaiting(int number, double probability) const
{
    int X = -1;
    double sum = 0;
    int iter = 0;
    do {
        sum += GeometricRand::variate(probability) + 1.0; // TODO: hash this one
        ++X;
    } while (sum <= number && ++iter < 1e9);
    return X;
}

double BinomialRand::variate() const
{
    if (npFloor <= generatorEdge)
    {
        if (p <= 0.5)
            return variateWaiting(n, p);
        return n - variateWaiting(n, q);
    }

    int Y = variateRejection();
    if (pRes > 0)
        Y += variateWaiting(n - Y, pRes / (1.0 - pFloor));
    return (p > 0.5) ? n - Y : Y;
}

std::complex<double> BinomialRand::CF(double t) const
{
    std::complex<double> y(q + p * std::cos(t), std::sin(t));
    return std::pow(y, n);
}

double BinomialRand::Median() const
{
    return std::round(np);
}

double BinomialRand::Mode() const
{
    return std::round(np + p - 1);
}

double BinomialRand::Skewness() const
{
    return (q - p) / std::sqrt(np * q);
}

double BinomialRand::ExcessKurtosis() const
{
    double y = 1.0 / (p * q);
    y -= 6.0;
    return y / n;
}
