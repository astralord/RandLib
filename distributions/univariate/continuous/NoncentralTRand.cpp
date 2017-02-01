#include "NoncentralTRand.h"
#include "NormalRand.h"

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
    logNu = std::log(nu);

    mu = noncentrality;
    cdfCoef = 0.5 * std::erfc(M_SQRT1_2 * mu); /// φ(-μ)

    T.SetDegree(nu);
}

double NoncentralTRand::Faux(double x, double nuAux, double muAux) const
{
    if (x == 0.0)
        return 0.0;
    double xSq = x * x;
    double temp = 0.5 * muAux * muAux;
    double logTemp = std::log(temp);
    double y = xSq / (xSq + nuAux);
    double halfNu = 0.5 * nuAux;
    double logY = std::log(y), aux = halfNu * std::log1p(-y);
    static constexpr double M_LN2_2 = 0.5 * M_LN2;

    /// go forward
    double sum1 = 0.0, sum2 = 0.0, addon1 = 0.0, addon2 = 0.0;
    double floorTemp = std::max(std::floor(temp), 1.0); /// starting point at most weight
    int j = floorTemp;
    double I1 = RandMath::regularizedBetaFun(y, j + 0.5, halfNu);
    double I2 = RandMath::regularizedBetaFun(y, j + 1.0, halfNu);
    double I1Temp = I1, I2Temp = I2;
    do {
        double jpHalf = j + 0.5, jp1 = j + 1.0;
        double p = j * logTemp - std::lgamma(jp1);
        p = std::exp(p);
        double q = j * logTemp - std::lgamma(j + 1.5) - M_LN2_2;
        q = std::exp(q);
        addon1 = p * I1;
        addon2 = q * I2;
        sum1 += addon1;
        sum2 += addon2;
        /// shift first regularized beta function
        double z = jpHalf * logY + aux;
        z = std::exp(z) / (jpHalf * RandMath::betaFun(jpHalf, halfNu));
        I1 -= z;
        /// shift second regularized beta function
        z = jp1 * logY + aux;
        z = std::exp(z) / (jp1 * RandMath::betaFun(jp1, halfNu));
        I2 -= z;
        ++j;
    } while (std::fabs(addon1) > MIN_POSITIVE * std::fabs(sum1) || std::fabs(addon2) > MIN_POSITIVE * std::fabs(sum2));

    double sum = sum1 + muAux * sum2;
    /// go backwards
    sum1 = 0.0, sum2 = 0.0;
    j = floorTemp - 1;
    I1 = I1Temp;
    I2 = I2Temp;
    addon1 = 1.0;
    addon2 = 1.0;
    while (j >= 0 && (std::fabs(addon1) > MIN_POSITIVE * std::fabs(sum1) || std::fabs(addon2) > MIN_POSITIVE * std::fabs(sum2))) {
        double jpHalf = j + 0.5, jp1 = j + 1.0;
        /// shift first regularized beta function
        double z = jpHalf * logY + aux;
        z = std::exp(z) / (jpHalf * RandMath::betaFun(jpHalf, halfNu));
        I1 += z;
        /// shift second regularized beta function
        z = jp1 * logY + aux;
        z = std::exp(z) / (jp1 * RandMath::betaFun(jp1, halfNu));
        I2 += z;
        double p = j * logTemp - std::lgamma(jp1);
        p = std::exp(p);
        double q = j * logTemp - std::lgamma(j + 1.5) - M_LN2_2;
        q = std::exp(q);
        addon1 = p * I1;
        addon2 = q * I2;
        sum1 += addon1;
        sum2 += addon2;
        --j;
    };
    sum += sum1 + muAux * sum2;
    sum *= 0.5;
    /// Avoid under- and overflow
    return (temp > 50 && sum > 0.0) ? std::exp(std::log(sum) - temp) : sum * std::exp(-temp);
}

double NoncentralTRand::f(double x) const
{
    if (mu == 0.0)
        return T.f(x);
    if (x == 0.0) {
        double y = std::lgamma(0.5 * nu + 0.5);
        double z = mu * mu + M_LNPI + logNu;
        y -= T.Y.GetLogGammaFunction();
        return std::exp(y - 0.5 * z);
    }
    double y = 0.0;
    if (x > 0.0)
        y = Faux(x * std::sqrt(1.0 + 2.0 / nu), nu + 2.0, mu) - Faux(x, nu, mu);
    else
        y = -Faux(x * std::sqrt(1.0 + 2.0 / nu), nu + 2.0, -mu) + Faux(x, nu, -mu);
    return nu * y / x;
}

double NoncentralTRand::F(double x) const
{
    if (mu == 0.0)
        return T.F(x);
    double y = ((x >= 0.0) ? Faux(x, nu, mu) : -Faux(x, nu, -mu));
    return cdfCoef + y;
}

double NoncentralTRand::Variate() const
{
    double X = NormalRand::StandardVariate() + mu;
    X /= std::sqrt(T.Y.Variate() / nu);
    // TODO: investigate if we can use Nakagami distribution instead
    return X;
}

double NoncentralTRand::Mean() const
{
    if (nu <= 1)
        return NAN;
    if (mu == 0.0)
        return 0.0;
    double mean = std::lgamma(0.5 * nu - 0.5);
    mean -= T.Y.GetLogGammaFunction();
    mean += 0.5 * (logNu - M_LN2);
    return mu * std::exp(mean);
}

double NoncentralTRand::Variance() const
{
    if (nu <= 2)
        return (nu > 1) ? INFINITY : NAN;
    double mean = Mean();
    return nu * (1.0 + mu * mu) / (nu - 2) + mean * mean;
}

double NoncentralTRand::Skewness() const
{
    if (nu <= 3)
        return NAN;
    double mean = Mean();
    double var = nu * (1.0 + mu * mu) / (nu - 2) + mean * mean;
    double thirdMoment = std::lgamma(0.5 * nu - 1.5);
    thirdMoment -= T.Y.GetLogGammaFunction();
    thirdMoment += 1.5 * logNu + 0.5 * M_LN2;
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
    thirdMoment += 1.5 * logNu + 0.5 * M_LN2;
    thirdMoment = 0.25 * std::exp(thirdMoment);
    thirdMoment *= mu * (3.0 + mu * mu);
    double mean = Mean();
    double kurtosis = fourthMoment - 4 * mean * thirdMoment;
    kurtosis += 6 * mean * mean * nu * (1.0 + mu * mu) / (nu - 2);
    kurtosis -= 3 * std::pow(mean, 4);
    double var = Variance();
    return kurtosis / (var * var);
}
