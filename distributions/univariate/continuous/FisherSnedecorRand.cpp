#include "FisherSnedecorRand.h"

FisherSnedecorRand::FisherSnedecorRand(int degree1, int degree2)
{
    SetDegrees(degree1, degree2);
}

std::string FisherSnedecorRand::Name() const
{
    return "Fisher(" + toStringWithPrecision(GetFirstDegree()) + ", " + toStringWithPrecision(GetSecondDegree()) + ")";
}

void FisherSnedecorRand::SetDegrees(int degree1, int degree2)
{
    d1 = std::max(degree1, 1);
    d2 = std::max(degree2, 1);

    B.SetParameters(0.5 * d1, 0.5 * d2);

    a = 0.5 * d1 - 1;
    d1_d2 = static_cast<double>(d1) / d2;
    c = -0.5 * (d1 + d2);
    d2_d1 = 1.0 / d1_d2;

    pdfCoef = (a + 1) * std::log(d1_d2);
    pdfCoef -= B.GetLogBetaFunction();
}

double FisherSnedecorRand::f(const double & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (a == 0.0)
            return std::exp(pdfCoef);
        return (a > 0) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

double FisherSnedecorRand::logf(const double & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (a == 0.0)
            return pdfCoef;
        return (a > 0) ? -INFINITY : INFINITY;
    }
    double y = a * std::log(x);
    y += c * std::log1p(d1_d2 * x);
    return pdfCoef + y;
}

double FisherSnedecorRand::F(const double & x) const
{
    return B.F(d1_d2 * x);
}

double FisherSnedecorRand::S(const double & x) const
{
    return B.S(d1_d2 * x);
}

double FisherSnedecorRand::Variate() const
{
    return d2_d1 * B.Variate();
}

void FisherSnedecorRand::Sample(std::vector<double> &outputData) const
{
    B.Sample(outputData);
    for (double &var : outputData)
        var = d2_d1 * var;
}

double FisherSnedecorRand::Mean() const
{
    return (d2 > 2) ? 1 + 2.0 / (d2 - 2) : INFINITY;
}

double FisherSnedecorRand::Variance() const
{
    if (d2 <= 4)
        return INFINITY;
    double variance = d2;
    variance /= d2 - 2;
    variance *= variance;
    variance *= 2 * (d1 + d2 - 2);
    variance /= d1;
    variance /= d2 - 4;
    return variance;
}

double FisherSnedecorRand::Median() const
{
    return d2_d1 * B.Median();
}

double FisherSnedecorRand::Mode() const
{
    if (d1 <= 2)
        return 0.0;
    return d2_d1 * (d1 - 2) / (d2 + 2);
}

double FisherSnedecorRand::Skewness() const
{
    if (d2 <= 6)
        return INFINITY;
    double skewness = 8.0 * (d2 - 4.0);
    double aux = d1 + d2 - 2;
    skewness /= d1 * aux;
    skewness = std::sqrt(skewness);
    skewness *= d1 + aux;
    skewness /= d2 - 6.0;
    return skewness;
}

double FisherSnedecorRand::ExcessKurtosis() const
{
    if (d2 <= 8)
        return INFINITY;
    double kurtosis = d2 - 2;
    kurtosis *= kurtosis;
    kurtosis *= d2 - 4;
    kurtosis /= d1;
    kurtosis /= d1 + d2 - 2;
    kurtosis += 5 * d2 - 22;
    kurtosis /= d2 - 6;
    kurtosis /= d2 - 8;
    return 12.0 * kurtosis;
}

double FisherSnedecorRand::quantileImpl(double p) const
{
  return d2_d1 * B.Quantile(p);
}

double FisherSnedecorRand::quantileImpl1m(double p) const
{
  return d2_d1 * B.Quantile1m(p);
}

std::complex<double> FisherSnedecorRand::CFImpl(double t) const
{
  return B.CF(d2_d1 * t);
}
