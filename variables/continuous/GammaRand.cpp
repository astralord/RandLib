#include "GammaRand.h"

GammaRand::GammaRand(double shape, double scale)
{
    setParameters(shape, scale);
}

std::string GammaRand::name()
{
    return "Gamma(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void GammaRand::setConstantsForGenerator()
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
    k = shape;
    if (k <= 0)
        k = 1.0;
    kInv = 1.0 / k;
    
    theta = scale;
    if (theta <= 0)
        theta = 1.0;
    thetaInv = 1.0 / theta;
   
    double k_round = std::round(k);
    if (RandMath::areEqual(k, k_round)) {
        k = k_round;
        cdfCoef = 1.0 / RandMath::factorial(k - 1);
    }
    else {
        cdfCoef = 1.0 / std::tgamma(k);
    }
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
    variateCoef = kInv + M_1_E;

    if (k > 3)
        setConstantsForGenerator();
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
    return (x < 0) ? 0 : cdfCoef * RandMath::lowerIncGamma(k, x * thetaInv);
}

double GammaRand::variate() const
{
    if (k < 5) {
        double k_round = std::round(k);
        if (RandMath::areEqual(k, k_round))
            return theta * variateForIntegerShape();
        if (RandMath::areEqual(k - 0.5, k_round))
            return theta * variateForIntegerShape();
        if (k <= 1)
            return theta * variateForSmallShape();
        if (k <= 3)
            return theta * variateForMediumShape();
    }
    return theta * variateForLargeShape();
}

void GammaRand::sample(QVector<double> &outputData) const
{
    if (k < 5) {
        double k_round = std::round(k);
        if (RandMath::areEqual(k, k_round)) {
            for (double &var : outputData)
                var = theta * variateForIntegerShape();
            return;
        }
        if (RandMath::areEqual(k - 0.5, k_round)) {
            for (double &var : outputData)
                var = theta * variateForHalfIntegerShape();
            return;
        }
        if (k <= 1) {
            for (double &var : outputData)
                var = theta * variateForSmallShape();
            return;
        }
        if (k <= 3) {
            for (double &var : outputData)
                var = theta * variateForMediumShape();
            return;
        }
    }
    
    for (double &var : outputData)
        var = theta * variateForLargeShape();
}

double GammaRand::Mean() const
{
    return k * theta;
}

double GammaRand::Variance() const
{
    return k * theta * theta;
}

double GammaRand::variateForIntegerShape() const
{
    double rv = 0;
    for (int i = 0; i < k; ++i)
        rv += ExponentialRand::standardVariate();
    return rv;
}

double GammaRand::variateForHalfIntegerShape() const
{
    double rv = 0;
    for (int i = 0; i < k - 1; ++i)
        rv += ExponentialRand::standardVariate();
    double N = NormalRand::standardVariate();
    return rv + .5 * N * N;
}

double GammaRand::variateForSmallShape() const
{
    double rv = 0;
    int iter = 0;
    do {
        double U = UniformRand::standardVariate();
        double p = k * variateCoef * U;
        double W = ExponentialRand::standardVariate();
        if (p <= 1)
        {
            rv = std::pow(p, kInv);
            if (rv <= W)
                return rv;
        }
        else
        {
            rv = -std::log(variateCoef * (1 - U));
            if ((1 - k) * std::log(rv) <= W)
                return rv;
        }
    } while (++iter < 1e9); /// one billion should be enough
    return NAN; /// shouldn't end up here
}

double GammaRand::variateForMediumShape() const
{
    double W1, W2;
    do {
        W1 = ExponentialRand::standardVariate();
        W2 = ExponentialRand::standardVariate();
    } while (W2 < (k - 1) * (W1 - std::log(W1) - 1));
    return k * W1;
}

double GammaRand::variateForLargeShape() const
{
    double rv = 0;
    int iter = 0;
    do {
        double U = UniformRand::standardVariate();
        if (U <= 0.0095722652)
        {
            double W1 = ExponentialRand::standardVariate();
            double W2 = ExponentialRand::standardVariate();
            rv = b * (1 + W1 / d);
            if (m * (rv / b - std::log(rv / m)) + c <= W2)
                return rv;
        }
        else
        {
            double N;
            do {
                N = NormalRand::standardVariate();
                rv = s * N + m;
            } while (rv < 0 || rv > b);
            U = UniformRand::standardVariate();
            double S = .5 * N * N;
            if (N > 0)
            {
                if (U < 1 - w * S)
                    return rv;
            }
            else if (U < 1 + S * (v * N - w))
                return rv;

            if (std::log(U) < m * std::log(rv / m) + m - rv + S)
                return rv;
        }
    } while (++iter < 1e9); /// one billion should be enough
    return NAN; /// shouldn't end up here
}

std::complex<double> GammaRand::CF(double t) const
{
    return std::pow(std::complex<double>(1.0, -theta * t), -k);
}

double GammaRand::Mode() const
{
    return (k < 1) ? 0 : (k - 1) * theta;
}

double GammaRand::Skewness() const
{
    return 2.0 / std::sqrt(k);
}

double GammaRand::ExcessKurtosis() const
{
    return 6.0 * kInv;
}

bool GammaRand::fitToData(const QVector<double> &sample)
{
    int N = sample.size();
    if (N == 0)
        return false;

    /// Calculate average
    long double average = 0.0L;
    long double logAverage = 0.0L;
    for (double var : sample) {
        if (var <= 0)
            return false;
        average += var;
        logAverage += std::log(var);
    }
    average /= N;
    logAverage /= N;

    /// Calculate initial guess for shape
    double s = std::log(average) - logAverage;
    double sm3 = s - 3.0, sp12 = 12.0 * s;
    double shape = sm3 * sm3 + sp12 + sp12;
    shape = std::sqrt(shape);
    shape -= sm3;
    shape /= sp12;

    if (!RandMath::findRoot([s] (double x)
    {
        return std::log(x) - RandMath::digamma(x) - s;
    },
    [] (double x)
    {
        return 1.0 / x - RandMath::trigamma(x);
    }, shape))
        return false;

    setParameters(shape, average / shape);
    return true;
}
