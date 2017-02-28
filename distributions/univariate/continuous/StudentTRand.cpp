#include "StudentTRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"

StudentTRand::StudentTRand(double degree, double location, double scale)
{
    SetDegree(degree);
    SetLocation(location);
    SetScale(scale);
}

std::string StudentTRand::Name() const
{
    if (mu == 0.0 && sigma == 1.0)
        return "Student's t(" + toStringWithPrecision(GetDegree()) + ")";
    return "Student's t(" + toStringWithPrecision(GetDegree()) + ", "
                          + toStringWithPrecision(GetLocation()) + ", "
                          + toStringWithPrecision(GetScale()) + ")";
}

void StudentTRand::SetDegree(double degree)
{
    nu = degree > 0 ? degree : 1;
    Y.SetParameters(0.5 * nu, 1.0);

    nup1Half = 0.5 * (nu + 1);
    pdfCoef = std::lgamma(nup1Half);
    pdfCoef -= 0.5 * M_LNPI;
    pdfCoef -= Y.GetLogGammaFunction();
    betaInv = std::exp(pdfCoef);
    pdfCoef -= 0.5 * std::log(nu);
}

void StudentTRand::SetLocation(double location)
{
    mu = location;
}

void StudentTRand::SetScale(double scale)
{
    sigma = scale > 0 ? scale : 1.0;
    logSigma = std::log(sigma);
}

double StudentTRand::f(double x) const
{
    /// adjustment
    double x0 = x - mu;
    x0 /= sigma;
    double xSq = x0 * x0;
    if (nu == 1) /// Cauchy distribution
        return M_1_PI / (sigma * (1 + xSq));
    if (nu == 2)
        return 1.0 / (sigma * std::pow(2.0 + xSq, 1.5));
    if (nu == 3) {
        double y = 3 + xSq;
        return 6 * M_SQRT3 * M_1_PI / (sigma * y * y);
    }
    double y = -nup1Half * std::log1p(xSq / nu);
    return std::exp(pdfCoef + y) / sigma;
}

double StudentTRand::logf(double x) const
{
    /// adjustment
    double x0 = x - mu;
    x0 /= sigma;
    double xSq = x0 * x0;
    if (nu == 1) /// Cauchy distribution
        return std::log(M_1_PI / (sigma * (1 + xSq)));
    if (nu == 2) {
        return -logSigma - 1.5 * std::log(2.0 + xSq);
    }
    if (nu == 3) {
        double y = 3 + xSq;
        y = 6 * M_SQRT3 * M_1_PI / (sigma * y * y);
        return std::log(y);
    }
    double y = -nup1Half * std::log1p(xSq / nu);
    return pdfCoef + y - logSigma;
}

double StudentTRand::F(double x) const
{
    double x0 = x - mu;
    x0 /= sigma;
    if (x0 == 0.0)
        return 0.5;
    if (nu == 1) {
        /// Cauchy distribution
        return 0.5 * M_1_PI * RandMath::atan(x);
    }
    if (nu == 2) {
        return 0.5 + 0.5 * x0 / std::sqrt(2 + x0 * x0);
    }
    if (nu == 3) {
        double y = M_SQRT3 * x0 / (x0 * x0 + 3);
        y += RandMath::atan(x0 / M_SQRT3);
        return 0.5 + M_1_PI * y;
    }
    double t = nu / (x0 * x0 + nu);
    double y = 0.5 * RandMath::incompleteBetaFun(t, 0.5 * nu, 0.5, 1.0 / betaInv) * betaInv;
    return (x0 > 0.0) ? (1 - y) : y;
}

double StudentTRand::S(double x) const
{
    double x0 = x - mu;
    x0 /= sigma;
    if (x0 == 0.0)
        return 0.5;
    if (nu == 1) {
        /// Cauchy distribution
        return 0.5 + M_1_PI * RandMath::atan(-x);
    }
    if (nu == 2) {
        return 0.5 - 0.5 * x0 / std::sqrt(2 + x0 * x0);
    }
    if (nu == 3) {
        double y = M_SQRT3 * x0 / (x0 * x0 + 3);
        y += RandMath::atan(x0 / M_SQRT3);
        return 0.5 - M_1_PI * y;
    }
    double t = nu / (x0 * x0 + nu);
    double y = 0.5 * RandMath::incompleteBetaFun(t, 0.5 * nu, 0.5, betaInv) * betaInv;
    return (x0 > 0.0) ? y : 1.0 - y;
}

double StudentTRand::Variate() const
{
    if (nu == 1)
        return CauchyRand::Variate(mu, sigma);
    return mu + sigma * NormalRand::StandardVariate() / Y.Variate();
}

void StudentTRand::Sample(std::vector<double> &outputData) const
{
    if (nu == 1) {
        for (double &var : outputData)
            var = CauchyRand::Variate(mu, sigma);
    }
    else {
        Y.Sample(outputData);
        for (double &var : outputData)
            var = mu + sigma * NormalRand::StandardVariate() / var;
    }
}

double StudentTRand::Mean() const
{
    return (nu > 1) ? mu : NAN;
}

double StudentTRand::Variance() const
{
    if (nu > 2)
        return sigma * sigma * nu / (nu - 2);
    return (nu > 1) ? INFINITY : NAN;
}

std::complex<double> StudentTRand::CFImpl(double t) const
{
    double x = std::sqrt(nu) * std::fabs(t * sigma); // value of sqrt(nu) can be hashed
    double vHalf = 0.5 * nu;
    double y = vHalf * std::log(x);
    y -= Y.GetLogGammaFunction();
    y -= (vHalf - 1) * M_LN2;
    y += RandMath::logModifiedBesselSecondKind(x, vHalf);
    double costmu = std::cos(t * mu), sintmu = std::sin(t * mu);
    std::complex<double> cf(costmu, sintmu);
    return std::exp(y) * cf;
}

double StudentTRand::quantileImpl(double p) const
{
    double temp = p - 0.5;
    if (nu == 1)
        return std::tan(M_PI * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (nu == 2)
        return sigma * 2.0 * temp * std::sqrt(0.5 / pq) + mu;
    if (nu == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * RandMath::sign(temp) * std::sqrt(beta);
    }
    // TODO: implement inverse of beta function
    return ContinuousDistribution::quantileImpl(p);
}

double StudentTRand::quantileImpl1m(double p) const
{
    double temp = 0.5 - p;
    if (nu == 1)
        return std::tan(M_PI * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (nu == 2)
        return sigma * 2 * temp * std::sqrt(0.5 / pq) + mu;
    if (nu == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * RandMath::sign(temp) * std::sqrt(beta);
    }
    return ContinuousDistribution::quantileImpl1m(p);
}

double StudentTRand::Median() const
{
    return mu;
}

double StudentTRand::Mode() const
{
    return mu;
}

double StudentTRand::Skewness() const
{
    return (nu > 3) ? 0.0 : NAN;
}

double StudentTRand::ExcessKurtosis() const
{
    if (nu > 4)
        return 6 / (nu - 4);
    return (nu > 2) ? INFINITY : NAN;
}
