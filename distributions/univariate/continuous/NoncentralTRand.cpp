#include "NoncentralTRand.h"
#include "NormalRand.h"
#include <functional>

NoncentralTRand::NoncentralTRand(double degree, double noncentrality) :
    T(degree)
{
    SetParameters(degree, noncentrality);
}

std::string NoncentralTRand::Name() const
{
    return "Noncentral-t(" + toStringWithPrecision(GetDegree()) + ", "
            + toStringWithPrecision(GetNoncentrality()) + ")";
}

void NoncentralTRand::SetParameters(double degree, double noncentrality)
{
    nu = (degree > 0.0) ? degree : 1.0;
    T.SetDegree(nu);
    sqrt1p2oNu = std::exp(0.5 * std::log1p(2.0 / nu));

    mu = noncentrality;
    PhiMu = 0.5 * std::erfc(-M_SQRT1_2 * mu); /// φ(μ)
    PhimMu = 0.5 * std::erfc(M_SQRT1_2 * mu); /// φ(-μ)

    /// precalculate values for pdf/cdf integration
    static constexpr double epsilon = 1e-16;
    ChiSquaredRand X(nu);
    nuCoefs.halfNu = 0.5 * nu;
    nuCoefs.qEpsCoef = std::sqrt(X.Quantile(epsilon) / nu);
    nuCoefs.q1mEpsCoef = std::sqrt(X.Quantile1m(epsilon) / nu);
    nuCoefs.logHalfNu = std::log(nuCoefs.halfNu);
    nuCoefs.lgammaHalfNu = T.Y.GetLogGammaFunction();
    double nup2 = nu + 2;
    X.SetDegree(nup2);
    nup2Coefs.halfNu = nuCoefs.halfNu + 1.0;
    nup2Coefs.qEpsCoef = std::sqrt(X.Quantile(epsilon) / nup2);
    nup2Coefs.q1mEpsCoef = std::sqrt(X.Quantile1m(epsilon) / nup2);
    nup2Coefs.logHalfNu = std::log1p(nuCoefs.halfNu);
    nup2Coefs.lgammaHalfNu = std::lgamma(nup2Coefs.halfNu);
}

DoublePair NoncentralTRand::getIntegrationLimits(double x, double muAux, const NoncentralTRand::nuStruct &nuAuxCoef) const
{
    static constexpr double normQuant = -37.5194;
    double A0 = std::max(-muAux, normQuant), B0 = -normQuant;
    double qEps = x * nuAuxCoef.qEpsCoef - muAux;
    double q1mEps = x * nuAuxCoef.q1mEpsCoef - muAux;
    return std::make_pair(std::max(qEps, A0), std::min(q1mEps, B0));
}

double NoncentralTRand::cdf(const double &x, const NoncentralTRand::nuStruct &nuAuxCoef, bool isCompl) const
{
    if (x < -mu)
        return upperTail(-x, -mu, nuAuxCoef, !isCompl);
    if (x < 0.0)
        return lowerTail(-x, -mu, nuAuxCoef, !isCompl);
    return (x < mu) ? lowerTail(x, mu, nuAuxCoef, isCompl) : upperTail(x, mu, nuAuxCoef, isCompl);
}

double NoncentralTRand::g(double z, const double &x, const nuStruct &nuAuxCoef, double muAux, bool lower) const
{
    if (z == -muAux)
        return 0.0;
    double y = -0.5 * z * z;
    y -= 0.5 * (M_LN2 + M_LNPI);
    double b = (z + muAux) / x;
    b *= b;
    b *= nuAuxCoef.halfNu;
    /// 'lower' stands for lower tail P(X < x),
    /// however in that case we need to call upper incomplete gamma function
    /// and vice versa
    y += lower ? RandMath::lqgamma(nuAuxCoef.halfNu, b, nuAuxCoef.logHalfNu, nuAuxCoef.lgammaHalfNu) :
                 RandMath::lpgamma(nuAuxCoef.halfNu, b, nuAuxCoef.logHalfNu, nuAuxCoef.lgammaHalfNu);
    return std::exp(y);
}

double NoncentralTRand::findMode(const double &x, double halfNuAux, double muAux, double A, double B) const
{
    /// find approximate value of mode
    double temp = (halfNuAux > 1) ? 8 * halfNuAux - 8 : 4.0;
    double mode = muAux * muAux + temp;
    double xSq = x * x;
    mode *= xSq;
    mode += 2 * temp * halfNuAux;
    mode = x * std::sqrt(mode);
    mode -= muAux * (xSq + 4 * halfNuAux);
    mode /= 2 * xSq + 4 * halfNuAux;
    /// sanity check
    if (mode <= A || mode >= B)
        mode = 0.5 * (A + B);
    return mode;
}

double NoncentralTRand::lowerTail(const double &x, double muAux, const nuStruct &nuAuxCoef, bool isCompl) const
{
    DoublePair intLimits = getIntegrationLimits(x, muAux, nuAuxCoef);
    double A = intLimits.first, B = intLimits.second;
    /// find peak of the integrand
    double mode = findMode(x, nuAuxCoef.halfNu, muAux, A, B);
    /// calculate two integrals
    std::function<double (double)> integrandPtr = std::bind(&NoncentralTRand::g, this, std::placeholders::_1, x, nuAuxCoef, muAux, true);
    double I1 = RandMath::integral(integrandPtr, A, mode);
    double I2 = RandMath::integral(integrandPtr, mode, B);
    double I = I1 + I2;
    /// for S(x) return 1 - F(x) w/o losing precision
    if (isCompl) {
        double PhimA = (A == -mu) ? PhiMu : 0.5 * std::erfc(M_SQRT1_2 * A); /// φ(-A)
        return PhimA - I;
    }
    double PhiA = (A == -mu) ? PhimMu : 0.5 * std::erfc(-M_SQRT1_2 * A); /// φ(A)
    return PhiA + I;
}

double NoncentralTRand::upperTail(const double &x, double muAux, const nuStruct &nuAuxCoef, bool isCompl) const
{
    DoublePair intLimits = getIntegrationLimits(x, muAux, nuAuxCoef);
    double A = intLimits.first, B = intLimits.second;
    /// find peak of the integrand
    double mode = findMode(x, nuAuxCoef.halfNu, muAux, A, B);
    /// calculate two integrals
    std::function<double (double)> integrandPtr = std::bind(&NoncentralTRand::g, this, std::placeholders::_1, x, nuAuxCoef, muAux, false);
    double I1 = RandMath::integral(integrandPtr, A, mode);
    double I2 = RandMath::integral(integrandPtr, mode, B);
    double I = I1 + I2;
    /// for S(x) return 1 - F(x) w/o losing precision
    if (isCompl) {
        double PhimB = (B == -mu) ? PhiMu : 0.5 * std::erfc(M_SQRT1_2 * B); /// φ(-B)
        return PhimB + I;
    }
    double PhiB = (B == -mu) ? PhimMu : 0.5 * std::erfc(-M_SQRT1_2 * B); /// φ(B)
    return PhiB - I;
}

double NoncentralTRand::f(const double & x) const
{
    if (mu == 0.0)
        return T.f(x);
    if (x == 0.0) {
        return std::exp(logf(x));
    }
    double y = cdf(x * sqrt1p2oNu, nup2Coefs, false);
    y -= cdf(x, nuCoefs, false);
    return nu * y / x;
}

double NoncentralTRand::logf(const double & x) const
{
    if (mu == 0.0)
        return T.logf(x);
    if (x == 0.0)
    {
        double y = mu * mu + M_LN2 + nuCoefs.logHalfNu;
        y *= -0.5;
        y -= T.logBetaFun;
        return y;
    }
    return std::log(f(x));
}

double NoncentralTRand::F(const double & x) const
{
    if (mu == 0.0)
        return T.F(x);
    if (x == 0.0)
        return PhimMu;
    return cdf(x, nuCoefs, false);
}

double NoncentralTRand::S(const double &x) const
{
    if (mu == 0.0)
        return T.S(x);
    if (x == 0.0)
        return PhiMu;
    return cdf(x, nuCoefs, true);
}

double NoncentralTRand::Variate() const
{
    double X = NormalRand::StandardVariate() + mu;
    X /= T.Y.Variate();
    return X;
}

void NoncentralTRand::Sample(std::vector<double> &outputData) const
{
    if (mu == 0.0)
        return T.Sample(outputData);
    T.Y.Sample(outputData);
    for (double &var : outputData)
        var = (mu + NormalRand::StandardVariate()) / var;
}

double NoncentralTRand::Mean() const
{
    if (nu <= 1)
        return NAN;
    if (mu == 0.0)
        return 0.0;
    double mean = std::lgamma(0.5 * nu - 0.5);
    mean -= T.Y.GetLogGammaFunction();
    mean += 0.5 * (nuCoefs.logHalfNu);
    return mu * std::exp(mean);
}

double NoncentralTRand::Variance() const
{
    if (nu <= 2)
        return (nu > 1) ? INFINITY : NAN;
    double mean = Mean();
    return nu * (1.0 + mu * mu) / (nu - 2) + mean * mean;
}

double NoncentralTRand::Mode() const
{
    if (mu == 0.0)
        return 0.0;
    double left = std::sqrt(nu / (nu + 2.5)); /// left boundary for mode / μ
    double right = std::sqrt(nu / (nu + 1.0)); /// right boundary for mode / μ
    double guess = 0.5 * mu * (left + right);
    double root = 0;
    RandMath::findMin([this] (double x)
    {
        return -f(x);
    }, guess, root);
    return root;
}

double NoncentralTRand::Skewness() const
{
    if (nu <= 3)
        return NAN;
    double mean = Mean();
    double var = nu * (1.0 + mu * mu) / (nu - 2) + mean * mean;
    double thirdMoment = std::lgamma(0.5 * nu - 1.5);
    thirdMoment -= T.Y.GetLogGammaFunction();
    thirdMoment += 1.5 * (nuCoefs.logHalfNu + M_LN2);
    thirdMoment = 0.25 * std::exp(thirdMoment);
    thirdMoment *= mu * (3.0 + mu * mu);
    double denominator = std::pow(var, 1.5);
    return (thirdMoment - 3 * mean * var - std::pow(mean, 3)) / denominator;
}

double NoncentralTRand::ExcessKurtosis() const
{
    if (nu <= 4)
        return (nu > 2) ? INFINITY : NAN;
    double fourthMoment = nu * nu / ((nu - 2) * (nu - 4));
    fourthMoment *= (std::pow(mu, 4) + 6 * mu * mu + 3);
    double thirdMoment = std::lgamma(0.5 * nu - 1.5);
    thirdMoment -= T.Y.GetLogGammaFunction();
    thirdMoment += 1.5 * (nuCoefs.logHalfNu + M_LN2);
    thirdMoment = 0.25 * std::exp(thirdMoment);
    thirdMoment *= mu * (3.0 + mu * mu);
    double mean = Mean();
    double kurtosis = fourthMoment - 4 * mean * thirdMoment;
    kurtosis += 6 * mean * mean * nu * (1.0 + mu * mu) / (nu - 2);
    kurtosis -= 3 * std::pow(mean, 4);
    double var = Variance();
    return kurtosis / (var * var);
}
