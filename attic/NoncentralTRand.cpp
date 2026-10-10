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

    halfMuSq = 0.5 * mu * mu;
    startingPoint = std::floor(halfMuSq);
    logHalfMuSq = 2 * std::log(std::fabs(mu)) - M_LN2;
    lgammaStartingPointpHalf = std::lgamma(startingPoint + 0.5);
    lgammaStartingPointp1 = RandMath::lfact(startingPoint);

    /// precalculate values for pdf/cdf calculation
    nuCoefs.halfNu = 0.5 * nu;
    nuCoefs.logHalfNu = std::log(nuCoefs.halfNu);
    nuCoefs.lgammaHalfNu = std::lgamma(nuCoefs.halfNu);
    nuCoefs.lgamma1 = std::lgamma(startingPoint + 0.5 + nuCoefs.halfNu);
    nuCoefs.lgamma2 = std::lgamma(startingPoint + 1 + nuCoefs.halfNu);
    nup2Coefs.halfNu = nuCoefs.halfNu + 1.0;
    nup2Coefs.logHalfNu = std::log1p(nuCoefs.halfNu);
    nup2Coefs.lgammaHalfNu = std::lgamma(nup2Coefs.halfNu);
    nup2Coefs.lgamma1 = std::lgamma(startingPoint + 0.5 + nup2Coefs.halfNu);
    nup2Coefs.lgamma2 = std::lgamma(startingPoint + 1 + nup2Coefs.halfNu);
}

double NoncentralTRand::cdfSeries(const double &x, const nuStruct &degreeCoef, double noncentrality) const
{
    if (x == 0.0)
        return 0.0;
    long double xSq = x * x;
    long double y = xSq / (xSq + 2 * degreeCoef.halfNu);
    long double logY = std::log(y), log1mY = std::log1p(-y);
    static constexpr double M_LN2_2 = 0.5 * M_LN2;

    /// go forward
    long double sum1 = 0.0, addon = 0.0;
    int j = startingPoint;
    long double logBetaFun1 = lgammaStartingPointpHalf + degreeCoef.lgammaHalfNu - degreeCoef.lgamma1;
    long double logBetaFun2 = lgammaStartingPointp1 + degreeCoef.lgammaHalfNu - degreeCoef.lgamma2;
    long double I1 = RandMath::ibeta(y, j + 0.5, degreeCoef.halfNu, logBetaFun1, logY, log1mY);
    long double I2 = RandMath::ibeta(y, j + 1.0, degreeCoef.halfNu, logBetaFun2, logY, log1mY);
    long double lgammajpHalf = lgammaStartingPointpHalf, lgammajp1 = lgammaStartingPointp1;
    long double lgammaj1 = degreeCoef.lgamma1, lgammaj2 = degreeCoef.lgamma2;
    long double I1Temp = I1, I2Temp = I2;
    size_t iter = 0;
    static constexpr size_t MAX_ITER = 1e5;
    do {
        double jpHalf = j + 0.5;
        int jp1 = j + 1;
        long double logJpHalf = std::log(jpHalf), logJp1 = std::log(jp1);
        long double temp = j * logHalfMuSq - halfMuSq;
        long double p = temp - lgammajp1;
        p = std::exp(p);
        long double q = temp - logJpHalf - lgammajpHalf - M_LN2_2;
        q = noncentrality * std::exp(q);
        addon = p * I1 + q * I2;
        sum1 += addon;
        if (std::fabs(addon) < MIN_POSITIVE * std::fabs(sum1))
            break;
        /// shift first regularized beta function
        long double z = jpHalf * logY + degreeCoef.halfNu * log1mY;
        long double logBeta = lgammajpHalf + degreeCoef.lgammaHalfNu - lgammaj1;
        z = std::exp(z - logBeta);
        z /= jpHalf;
        I1 -= z;
        /// shift second regularized beta function
        z = jp1 * logY + degreeCoef.halfNu * log1mY;
        logBeta = lgammajp1 + degreeCoef.lgammaHalfNu - lgammaj2;
        z = std::exp(z - logBeta);
        z /= jp1;
        I2 -= z;
        /// go to the next term
        lgammajpHalf += logJpHalf;
        lgammajp1 += logJp1;
        lgammaj1 += std::log(jpHalf + degreeCoef.halfNu);
        lgammaj2 += std::log(jp1 + degreeCoef.halfNu);
        ++j;
    } while (++iter < MAX_ITER);

    /// go backwards
    long double sum2 = 0.0;
    j = startingPoint;
    I1 = I1Temp;
    I2 = I2Temp;
    lgammajpHalf = lgammaStartingPointpHalf;
    long double lgammajmHalf = lgammaStartingPointpHalf;
    long double lgammaj = lgammaStartingPointp1;
    lgammaj1 = degreeCoef.lgamma1;
    lgammaj2 = degreeCoef.lgamma2;
    while (j > 0) {
        /// move to j - 1
        double jmHalf = j - 0.5;
        lgammajmHalf -= std::log(jmHalf);
        lgammaj -= std::log(j);
        lgammaj1 -= std::log(jmHalf + degreeCoef.halfNu);
        lgammaj2 -= std::log(j + degreeCoef.halfNu);
        /// shift first regularized beta function
        long double z = jmHalf * logY + degreeCoef.halfNu * log1mY;
        long double logBeta = lgammajmHalf + degreeCoef.lgammaHalfNu - lgammaj1;
        z = std::exp(z - logBeta);
        z /= jmHalf;
        I1 += z;
        /// shift second regularized beta function
        z = j * logY + degreeCoef.halfNu * log1mY;
        logBeta = lgammaj + degreeCoef.lgammaHalfNu - lgammaj2;
        z = std::exp(z - logBeta);
        z /= j;
        I2 += z;
        long double temp = (j - 1) * logHalfMuSq - halfMuSq;
        long double p = temp - lgammaj; // p and q get by relation with previous values
        p = std::exp(p);
        long double q = temp - lgammajpHalf - M_LN2_2;
        q = noncentrality * std::exp(q);
        addon = p * I1 + q * I2;
        sum2 += addon;
        if (std::fabs(addon) < MIN_POSITIVE * std::fabs(sum2))
            break;
        --j;
        lgammajpHalf = lgammajmHalf;
    };
    //if (std::fabs(x) < 0.05)
        //qDebug() << "SUM:" << x << (double)sum1 << (double)sum2 << (double)(0.5 * (sum1 + sum2));
    return 0.5 * (sum1 + sum2);
}

double NoncentralTRand::cdfComplSeries(const double &x, const NoncentralTRand::nuStruct &degreeCoef, double noncentrality) const
{
    // TODO: find bug here, sum is always around 1e-15
    long double y = degreeCoef.halfNu / (degreeCoef.halfNu + 0.5 * x * x);
    double logY = std::log(y), log1mY = std::log1p(-y);
    static constexpr long double M_LN2_2 = 0.5 * M_LN2;

    /// go forward
    long double sum1 = 0.0, addon = 0.0;
    int j = startingPoint;
    long double logBetaFun1 = lgammaStartingPointpHalf + degreeCoef.lgammaHalfNu - degreeCoef.lgamma1;
    long double logBetaFun2 = lgammaStartingPointp1 + degreeCoef.lgammaHalfNu - degreeCoef.lgamma2;
    long double I1 = RandMath::ibeta(y, degreeCoef.halfNu, j + 0.5, logBetaFun1, logY, log1mY);
    long double I2 = RandMath::ibeta(y, degreeCoef.halfNu, j + 1.0, logBetaFun2, logY, log1mY);
    long double lgammajpHalf = lgammaStartingPointpHalf, lgammajp1 = lgammaStartingPointp1;
    long double lgammaj1 = degreeCoef.lgamma1, lgammaj2 = degreeCoef.lgamma2;
    long double I1Temp = I1, I2Temp = I2;
    size_t iter = 0;
    static constexpr size_t MAX_ITER = 1e5;
    do {
        double jpHalf = j + 0.5;
        int jp1 = j + 1;
        long double logJpHalf = std::log(jpHalf), logJp1 = std::log(jp1);
        long double temp = j * logHalfMuSq - halfMuSq;
        long double p = temp - lgammajp1;
        p = std::exp(p);
        long double q = temp - logJpHalf - lgammajpHalf - M_LN2_2;
        q = noncentrality * std::exp(q);
        addon = p * I1 + q * I2;
        sum1 += addon;
        if (std::fabs(addon) < MIN_POSITIVE * std::fabs(sum1))
            break;
        /// shift first regularized beta function
        long double z = degreeCoef.halfNu  * logY + jpHalf * log1mY;
        long double logBeta = lgammajpHalf + degreeCoef.lgammaHalfNu - lgammaj1;
        z = std::exp(z - logBeta);
        z /= jpHalf;
        I1 += z;
        /// shift second regularized beta function
        z = degreeCoef.halfNu * logY + jp1 * log1mY;
        logBeta = lgammajp1 + degreeCoef.lgammaHalfNu - lgammaj2;
        z = std::exp(z - logBeta);
        z /= jp1;
        I2 += z;
        /// go to the next term
        lgammajpHalf += logJpHalf;
        lgammajp1 += logJp1;
        lgammaj1 += std::log(jpHalf + degreeCoef.halfNu);
        lgammaj2 += std::log(jp1 + degreeCoef.halfNu);
        ++j;
    } while (++iter < MAX_ITER);

    /// go backwards
    long double sum2 = 0.0;
    j = startingPoint;
    I1 = I1Temp;
    I2 = I2Temp;
    lgammajpHalf = lgammaStartingPointpHalf;
    long double lgammajmHalf = lgammaStartingPointpHalf;
    long double lgammaj = lgammaStartingPointp1;
    lgammaj1 = degreeCoef.lgamma1;
    lgammaj2 = degreeCoef.lgamma2;
    while (j > 0) {
        /// move to j - 1
        double jmHalf = j - 0.5;
        lgammajmHalf -= std::log(jmHalf);
        lgammaj -= std::log(j);
        lgammaj1 -= std::log(jmHalf + degreeCoef.halfNu);
        lgammaj2 -= std::log(j + degreeCoef.halfNu);
        /// shift first regularized beta function
        long double z = degreeCoef.halfNu * logY + jmHalf * log1mY;
        long double logBeta = lgammajmHalf + degreeCoef.lgammaHalfNu - lgammaj1;
        z = std::exp(z - logBeta);
        z /= jmHalf;
        I1 -= z;
        /// shift second regularized beta function
        z = degreeCoef.halfNu * logY + j * log1mY;
        logBeta = lgammaj + degreeCoef.lgammaHalfNu - lgammaj2;
        z = std::exp(z - logBeta);
        z /= j;
        I2 -= z;
        long double temp = (j - 1) * logHalfMuSq - halfMuSq;
        long double p = temp - lgammaj;
        p = std::exp(p);
        long double q = temp - lgammajpHalf - M_LN2_2;
        q = noncentrality * std::exp(q);
        addon = p * I1 + q * I2;
        sum2 += addon;
        if (std::fabs(addon) < MIN_POSITIVE * std::fabs(sum2))
            break;
        --j;
        lgammajpHalf = lgammajmHalf;
    };
    //if (std::fabs(x) < 0.05)
        //qDebug() << "SUM:" << x << (double)sum1 << (double)sum2 << (double)(0.5 * (sum1 + sum2));
    return 0.5 * (sum1 + sum2);
}

double NoncentralTRand::logPdfAtZero() const
{
    double y = mu * mu + M_LN2 + nuCoefs.logHalfNu;
    y *= -0.5;
    y -= T.logBetaFun;
    return y;
}

double NoncentralTRand::pdfCommon(const double &x, double noncentrality) const
{
    double y = 0.0;
    if (x >= 0.0) {
        y = cdfSeries(x * sqrt1p2oNu, nup2Coefs, noncentrality);
        y -= cdfSeries(x, nuCoefs, noncentrality);
    }
    else {
        y = cdfComplSeries(-x * sqrt1p2oNu, nup2Coefs, -noncentrality);
        y -= cdfComplSeries(-x, nuCoefs, -noncentrality);
    }
    y *= nu / x;
    return std::max(y, 0.0);
}

double NoncentralTRand::f(const double & x) const
{
    if (mu == 0.0)
        return T.f(x);
    return (x == 0.0) ? std::exp(logPdfAtZero()) : pdfCommon(x, mu);
}

double NoncentralTRand::logf(const double & x) const
{
    if (mu == 0.0)
        return T.logf(x);
    if (x == 0.0)
        return logPdfAtZero();
    double y = pdfCommon(x, mu);
    return (y > 0.0) ? std::log(y) : -INFINITY;
}

double NoncentralTRand::F(const double & x) const
{
    if (mu == 0.0)
        return T.F(x);
    if (x == 0.0)
        return PhimMu;

    // mu > 0
    if (x >= 0.0)
        return PhimMu + cdfSeries(x, nuCoefs, mu);
    return PhimMu - cdfSeries(-x, nuCoefs, -mu);

    /*double y = (x >= 0.0) ? PhimMu + cdfSeries(x, nuCoefs, mu) : cdfComplSeries(-x, nuCoefs, -mu);
    return std::max(y, 0.0);*/
}

double NoncentralTRand::S(const double &x) const
{
    if (mu == 0.0)
        return T.S(x);
    double y = (x >= 0.0) ? cdfComplSeries(x, nuCoefs, mu) : PhiMu + cdfSeries(-x, nuCoefs, -mu);
    return std::max(y, 0.0);
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
    double temp = T.Y.GetLogGammaShapeRatio();
    temp += 0.5 * nuCoefs.logHalfNu;
    return 2 * mu / (nu - 1) * std::exp(temp);
}

double NoncentralTRand::Variance() const
{
    if (nu <= 2)
        return (nu > 1) ? INFINITY : NAN;
    double mean = Mean();
    return nu * (1.0 + mu * mu) / (nu - 2) + mean * mean; // TODO: Wikipedia lies, re-check
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
    double thirdMoment = std::lgamma(0.5 * nu - 1.5); //TODO: use hashed GetLogGammaShapeRatio
    thirdMoment -= T.Y.GetLogGammaFunction();
    thirdMoment += 1.5 * (nuCoefs.logHalfNu + M_LN2);
    thirdMoment = 0.25 * std::exp(thirdMoment);
    thirdMoment *= mu * (3.0 + mu * mu);
    double denominator = std::pow(var, 1.5);
    return (thirdMoment - 3 * mean * var - std::pow(mean, 3)) / denominator;  // TODO: Wikipedia lies, re-check
}

double NoncentralTRand::ExcessKurtosis() const
{
    if (nu <= 4)
        return (nu > 2) ? INFINITY : NAN;
    double fourthMoment = nu * nu / ((nu - 2) * (nu - 4));
    fourthMoment *= (std::pow(mu, 4) + 6 * mu * mu + 3);
    double thirdMoment = std::lgamma(0.5 * nu - 1.5); //TODO: use hashed GetLogGammaShapeRatio
    thirdMoment -= T.Y.GetLogGammaFunction();
    thirdMoment += 1.5 * (nuCoefs.logHalfNu + M_LN2);
    thirdMoment = 0.25 * std::exp(thirdMoment);
    thirdMoment *= mu * (3.0 + mu * mu);
    double mean = Mean();
    double kurtosis = fourthMoment - 4 * mean * thirdMoment;
    kurtosis += 6 * mean * mean * nu * (1.0 + mu * mu) / (nu - 2);
    kurtosis -= 3 * std::pow(mean, 4);
    double var = Variance();
    return kurtosis / (var * var);  // TODO: Wikipedia lies, re-check
}

double NoncentralTRand::quantileImpl(double p) const
{
    if (p > 1e-5)
        return ContinuousDistribution::quantileImpl(p);
    /// we need to be cautious with small values of p
    /// initial guess is obtained by inverse error function approximation
    double guess = -2 * std::erfc(2 * p);
    guess += mu;
    guess /= 1 - 0.25 * nu;
    double logP = std::log(p);
    if (RandMath::findRoot([this, logP] (double x)
    {
        double logCdf = std::log(F(x)), logPdf = logf(x);
        double first = logCdf - logP;
        double second = std::exp(logPdf - logCdf);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double NoncentralTRand::quantileImpl1m(double p) const
{
    if (p > 1e-5)
        return ContinuousDistribution::quantileImpl1m(p);
    /// we need to be cautious with small values of p
    /// initial guess is obtained by inverse error function approximation
    double guess = 2 * std::erfc(2 * p);
    guess += mu;
    guess /= 1 - 0.25 * nu;
    double logP = std::log(p);
    if (RandMath::findRoot([this, logP] (double x)
    {
        double logCcdf = std::log(S(x)), logPdf = logf(x);
        double first = logP - logCcdf;
        double second = std::exp(logPdf - logCcdf);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}
