#include "GammaRand.h"

GammaRand::GammaRand(double shape, double scale)
{
    setParameters(shape, scale);
}

void GammaRand::setConstants()
{
    m = k - 1;
    s_2 = std::sqrt(8.0 * k / 3) + k;
    s = std::sqrt(s_2);
    d = M_SQRT2 * M_SQRT3 * s_2;
    b = d + m;
    w = s_2 / (m - 1);
    v = (s_2 + s_2) / (m * std::sqrt(k));
    c = b + std::log(s * d / b) - m - m - 3.7203285;
}

void GammaRand::setParameters(double shape, double scale)
{
    k = std::max(shape, MIN_POSITIVE);
    kInv = 1.0 / k;
    theta = std::max(scale, MIN_POSITIVE);
    thetaInv = 1.0 / theta;

    cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
    valueCoef = kInv + M_1_E;

    if (k <= 1)
        U.setBoundaries(0, 1 + k * M_1_E);
    else if (k > 3)
        setConstants();
}

void GammaRand::setShape(double shape)
{
    k = std::max(shape, MIN_POSITIVE);
    kInv = 1.0 / k;

    cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
    valueCoef = kInv + M_1_E;

    if (k <= 1)
        U.setBoundaries(0, 1 + k * M_1_E);
    else if (k > 3)
        setConstants();
}

void GammaRand::setScale(double scale)
{
    theta = std::max(scale, MIN_POSITIVE);
    thetaInv = 1.0 / theta;
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
}

double GammaRand::f(double x) const
{
    if (x < 0)
        return 0;
    double y = std::pow(x, k - 1);
    y *= std::exp(-x * thetaInv);
    return pdfCoef * y;
}

double GammaRand::F(double x) const
{
    if (x < 0)
        return 0;
    return cdfCoef * RandMath::lowerIncGamma(k, x * thetaInv);
}

double GammaRand::value()
{
    double rv = 0;

    /// GA algorithm for integers and half-integers
    int k_int = (int)std::floor(k);
    if (std::fabs(k - k_int) < MIN_POSITIVE) /// Erlang distribution
    {
        for (int i = 0; i != k_int; ++i)
            rv += E.value();
        return theta * rv;
    }

    if (std::fabs(k - k_int - .5) < MIN_POSITIVE)
    {
        for (int i = 0; i != k_int; ++i)
            rv += E.value();
        double n = N.value();
        rv += .5 * n * n;
        return theta * rv;
    }

    /// GS algorithm for small k < 1
    if (k <= 1.0)
    {
        int iter = 0;
        do {
            double P = U.value();
            double e = E.value();
            if (P <= 1)
            {
                rv = std::pow(P, kInv);
                if (rv <= e)
                    return theta * rv;
            }
            else
            {
                rv = -std::log(valueCoef - kInv * P);
                if ((1 - k) * std::log(rv) <= e)
                    return theta * rv;
            }
        } while (++iter < 1e9); /// one billion should be enough
    }

    /// GP algorithm for 1 < k < 3
    if (k <= 3.0)
    {
        double e1, e2;
        do {
            e1 = E.value();
            e2 = E.value();
        } while (e2 < (k - 1) * (e1 - std::log(e1) - 1));
        return theta * k * e1;
    }

    /// GO algorithm for most common case k > 3
    int iter = 0;
    do {
        double u = U.value();
        if (u <= 0.0095722652)
        {
            double e1 = E.value();
            double e2 = E.value();
            rv = b * (1 + e1 / d);
            if (m * (rv / b - std::log(rv / m)) + c <= e2)
                return theta * rv;
        }
        else
        {
            double n;
            do {
                n = N.value();
                rv = s * n + m;
            } while (rv < 0 || rv > b);
            u = U.value();
            double S = .5 * n * n;
            if (n > 0)
            {
                if (u < 1 - w * S)
                    return theta * rv;
            }
            else if (u < 1 + S * (v * n - w))
                return theta * rv;

            if (std::log(u) < m * std::log(rv / m) + m - rv + S)
                return theta * rv;

        }
    } while (++iter < 1e9); /// one billion should be enough
    return 0;
}
