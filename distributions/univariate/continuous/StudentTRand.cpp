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
    v = degree > 0 ? degree : 1;
    Y.SetParameters(0.5 * v, 0.5);

    vp1Half = 0.5 * (v + 1);
    pdfCoef = std::lgamma(vp1Half);
    pdfCoef -= 0.5 * M_LNPI;
    pdfCoef -= Y.GetLogGammaFunction();
    betaInv = std::exp(pdfCoef);
    pdfCoef -= 0.5 * std::log(v);
}

void StudentTRand::SetLocation(double location)
{
    mu = location;
}

void StudentTRand::SetScale(double scale)
{
    sigma = scale > 0 ? scale : 1.0;
}

double StudentTRand::f(double x) const
{
    // TODO: consider special cases
    /// adjustment
    x -= mu;
    x /= sigma;

    double y = 1 + x * x / v;
    y = -vp1Half * std::log(y);
    return std::exp(pdfCoef + y) / sigma;
}

double StudentTRand::F(double x) const
{
    x -= mu;
    x /= sigma;
    if (x == 0.0)
        return 0.5;

    // TODO: consider v == 3 and v == 4
    if (v == 1)
        return 0.5 + std::atan(x) / M_PI;
    if (v == 2)
        return 0.5 * (1.0 + x / std::sqrt(2 + x * x));

    double t = v / (x * x + v);
    double y = 0.5 * RandMath::incompleteBetaFun(t, 0.5 * v, 0.5) * betaInv;
    return (x > 0.0) ? (1 - y) : y;
}

double StudentTRand::Variate() const
{
    //v = inf -> normal
    if (v == 1)
        return CauchyRand::Variate(mu, sigma);
    return mu + sigma * NormalRand::StandardVariate() / std::sqrt(Y.Variate() / v);
}

void StudentTRand::Sample(std::vector<double> &outputData) const
{
    //v = inf -> normal
    if (v == 1) {
        for (double &var : outputData)
            var = CauchyRand::Variate(mu, sigma);
    }
    else {
        Y.Sample(outputData);
        for (double &var : outputData)
            var = mu + sigma * NormalRand::StandardVariate() / std::sqrt(var / v);
    }
}

double StudentTRand::Mean() const
{
    return (v > 1) ? mu : NAN;
}

double StudentTRand::Variance() const
{
    if (v > 2)
        return sigma * sigma * v / (v - 2);
    return (v > 1) ? INFINITY : NAN;
}

std::complex<double> StudentTRand::CF(double t) const
{
    double x = std::sqrt(v) * std::fabs(t * sigma); // value of sqrt(v) can be hashed
    double vHalf = 0.5 * v;
    double y = vHalf * std::log(x);
    y -= Y.GetLogGammaFunction();
    y -= (vHalf - 1) * M_LN2;
    std::complex<double> cf(y, t * mu);
    cf = std::exp(cf);
    cf *= RandMath::modifiedBesselSecondKind(x, vHalf);
    return cf;
}

double StudentTRand::quantileImpl(double p) const
{
    // TODO: consider v == 3
    double temp = p - 0.5;
    if (v == 1)
        return std::tan(M_PI * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (v == 2)
        return sigma * 2.0 * temp * std::sqrt(0.5 / pq) + mu;
    if (v == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * RandMath::sign(temp) * std::sqrt(beta);
    }
    return ContinuousDistribution::quantileImpl(p);
}

double StudentTRand::quantileImpl1m(double p) const
{
    // TODO: consider v == 3
    double temp = 0.5 - p;
    if (v == 1)
        return std::tan(M_PI * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (v == 2)
        return sigma * 2 * temp * std::sqrt(0.5 / pq) + mu;
    if (v == 4)
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
    return (v > 3) ? 0.0 : NAN;
}

double StudentTRand::ExcessKurtosis() const
{
    if (v > 4)
        return 6 / (v - 4);
    return (v > 2) ? INFINITY : NAN;
}
