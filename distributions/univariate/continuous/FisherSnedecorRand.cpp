#include "FisherSnedecorRand.h"

FisherSnedecorRand::FisherSnedecorRand(int degree1, int degree2)
{
    setDegrees(degree1, degree2);
}

std::string FisherSnedecorRand::name() const
{
    return "Fisher(" + toStringWithPrecision(getFirstDegree()) + ", " + toStringWithPrecision(getSecondDegree()) + ")";
}

void FisherSnedecorRand::setDegrees(int degree1, int degree2)
{
    d1 = std::max(degree1, 1);
    d2 = std::max(degree2, 1);

    B.setParameters(.5 * d1, .5 * d2);

    a = .5 * d1 - 1;
    d1_d2 = static_cast<double>(d1) / d2;
    c = -.5 * (d1 + d2);
    d2_d1 = 1.0 / d1_d2;

    pdfCoef = (a + 1) * std::log(d1_d2);
    pdfCoef -= B.getLogBetaFunction();
}

double FisherSnedecorRand::f(double x) const
{
    if (x < 0)
        return 0.0;
    if (x == 0) {
        if (a == 0)
            return std::exp(pdfCoef);
        return (a > 0) ? 0.0 : INFINITY;
    }
    double y = a * std::log(x);
    y += c * std::log1p(d1_d2 * x);
    return std::exp(pdfCoef + y);
}

double FisherSnedecorRand::F(double x) const
{
    return B.F(d1_d2 * x);
}

double FisherSnedecorRand::variate() const
{
    return d2_d1 * B.variate();
}

std::complex<double> FisherSnedecorRand::CF(double t) const
{
    return B.CF(d2_d1 * t);
}

void FisherSnedecorRand::sample(std::vector<double> &outputData) const
{
    B.sample(outputData);
    for (double &var : outputData)
        var = d2_d1 * var;
}

double FisherSnedecorRand::Mean() const
{
    return (d2 > 2) ? d2 / (d2 - 2) : INFINITY;
}

double FisherSnedecorRand::Variance() const
{
    if (d2 <= 4)
        return INFINITY;
    double numerator = 2 * d2 * d2 * (d1 + d2 - 2);
    double denominator = d2 - 2;
    denominator *= denominator;
    denominator *= d1 * (d2 - 4);
    return numerator / denominator;
}

double FisherSnedecorRand::QuantileImpl(double p) const
{
    return d2_d1 * B.QuantileImpl(p);
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
    double skewness = 8 * (d2 - 4);
    double aux = d1 + d2 - 2;
    skewness /= d1 * aux;
    skewness = std::sqrt(skewness);
    skewness *= d1 + aux;
    skewness /= d2 - 6;
    return skewness;
}

double FisherSnedecorRand::ExcessKurtosis() const
{
    if (d2 <= 8)
        return INFINITY;
    double kurtosis = d2 - 2;
    kurtosis *= kurtosis;
    kurtosis *= d2 - 4;
    kurtosis /= d1 * (d1 + d2 - 2);
    kurtosis += 5 * d2 - 22;
    kurtosis /= (d2 - 6) * (d2 - 8);
    return 12.0 * kurtosis;
}
