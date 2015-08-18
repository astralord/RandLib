#include "GammaRand.h"

GammaRand::GammaRand(double shape, double scale)
{
    setParameters(shape, scale);
}

void GammaRand::setName()
{
    nameStr = "Gamma(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
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

    if (std::fabs(k - std::round(k)) < MIN_POSITIVE)
        cdfCoef = 1.0 / RandMath::fastFactorial(k);
    else
        cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
    valueCoef = kInv + M_1_E;

    if (k <= 1) {
        U.setBoundaries(0, 1 + k * M_1_E);
    }
    else if (k > 3)
    {
        setConstants();
        U.setBoundaries(0, 1);
    }
    setName();
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

double GammaRand::variate() const
{
    double rv = 0;
    int k_int = static_cast<int>(k);
    if (std::fabs(k - k_int) < MIN_POSITIVE) {
        rv = variateForIntegerShape();
    }
    else if (std::fabs(k - k_int - .5) < MIN_POSITIVE) {
        rv = variateForHalfIntegerShape();
    }
    else if (k <= 1) {
        rv = variateForSmallShape();
    }
    else if (k <= 3) {
        rv = variateForMediumShape();
    }
    else {
        rv = variateForLargeShape();
    }

    return theta * rv;
}

void GammaRand::sample(QVector<double> &outputData)
{
    int k_int = static_cast<int>(k);

    if (std::fabs(k - k_int) < MIN_POSITIVE) {
        for (double &var : outputData)
            var = theta * variateForIntegerShape();
    }
    else if (std::fabs(k - k_int - .5) < MIN_POSITIVE) {
        for (double &var : outputData)
            var = theta * variateForHalfIntegerShape();
    }
    else if (k <= 1) {
        for (double &var : outputData)
            var = theta * variateForSmallShape();
    }
    else if (k <= 3) {
        for (double &var : outputData)
            var = theta * variateForMediumShape();
    }
    else {
        for (double &var : outputData)
            var = theta * variateForLargeShape();
    }
}

double GammaRand::variateForIntegerShape() const
{
    double rv = 0;
    for (int i = 0; i < k; ++i)
        rv += W.variate();
    return rv;
}

double GammaRand::variateForHalfIntegerShape() const
{
    double rv = 0;
    for (int i = 0; i < k - 1; ++i)
        rv += W.variate();
    double n = N.variate();
    return rv + .5 * n * n;
}

double GammaRand::variateForSmallShape() const
{
    double rv = 0;
    int iter = 0;
    do {
        double P = U.variate();
        double e = W.variate();
        if (P <= 1)
        {
            rv = std::pow(P, kInv);
            if (rv <= e)
                return rv;
        }
        else
        {
            rv = -std::log(valueCoef - kInv * P);
            if ((1 - k) * std::log(rv) <= e)
                return rv;
        }
    } while (++iter < 1e9); /// one billion should be enough
    return 0; /// shouldn't end up here
}

double GammaRand::variateForMediumShape() const
{
    double E1, E2;
    do {
        E1 = W.variate();
        E2 = W.variate();
    } while (E2 < (k - 1) * (E1 - std::log(E1) - 1));
    return k * E1;
}

double GammaRand::variateForLargeShape() const
{
    double rv = 0;
    int iter = 0;
    do {
        double u = U.variate();
        if (u <= 0.0095722652)
        {
            double e1 = W.variate();
            double e2 = W.variate();
            rv = b * (1 + e1 / d);
            if (m * (rv / b - std::log(rv / m)) + c <= e2)
                return rv;
        }
        else
        {
            double n;
            do {
                n = N.variate();
                rv = s * n + m;
            } while (rv < 0 || rv > b);
            u = U.variate();
            double S = .5 * n * n;
            if (n > 0)
            {
                if (u < 1 - w * S)
                    return rv;
            }
            else if (u < 1 + S * (v * n - w))
                return rv;

            if (std::log(u) < m * std::log(rv / m) + m - rv + S)
                return rv;
        }
    } while (++iter < 1e9); /// one billion should be enough
    return 0; /// shouldn't end up here
}
