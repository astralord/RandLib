#include "RayleighRand.h"

RayleighRand::RayleighRand(double scale)
{
    setScale(scale);
}

void RayleighRand::setScale(double scale)
{
    sigma = std::max(scale, MIN_POSITIVE);
    sigmaSqInv_2 = 0.5 / sigma * sigma;
}

double RayleighRand::f(double x) const
{
    if (x <= 0)
        return 0.0;
    double y = x * sigmaSqInv_2;
    y *= std::exp(-x * x * sigmaSqInv_2);
    return y + y;
}

double RayleighRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    return 1.0 - std::exp(-x * x * sigmaSqInv_2);
}

double RayleighRand::variate()
{
    return std::sqrt(ExponentialRand::variate(sigmaSqInv_2));
}

bool RayleighRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate sigma^2
    double sigmaEst = 0.0;
    int N = sample.size();
    for (double var : sample)
    {
        if (var < 0)
            return false;
        sigmaEst += var * var;
    }
    sigmaEst *= .5 / N;

    /// Calculate unbiased sigma
    sigmaEst = std::sqrt(sigmaEst);

    if (N > 30)
        setScale((1 + 0.1252 / N) * sigmaEst); /// err < 1e-6
    else
    {
        double coef = RandMath::fastFactorial(N - 1);
        coef *= N * coef;
        coef *= M_1_SQRTPI * std::sqrt(static_cast<double>(N));
        coef /= RandMath::fastFactorial(N << 1);
        int pow2N = 1 << N; /// < 2^31
        coef *= pow2N;
        coef *= pow2N;

        setScale(coef * sigmaEst);
    }
    return true;
}
