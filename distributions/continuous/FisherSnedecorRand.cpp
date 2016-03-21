#include "FisherSnedecorRand.h"

FisherSnedecorRand::FisherSnedecorRand(int degree1, int degree2)
{
    setDegrees(degree1, degree2);
}

std::string FisherSnedecorRand::name()
{
    return "Fisher(" + toStringWithPrecision(getFirstDegree()) + ", " + toStringWithPrecision(getSecondDegree()) + ")";
}

void FisherSnedecorRand::setDegrees(int degree1, int degree2)
{
    d1 = std::max(degree1, 1);
    d2 = std::max(degree2, 1);

    B.setShapes(.5 * d1, .5 * d2);

    a = .5 * d1 - 1;
    d1_d2 = static_cast<double>(d1) / d2;
    c = -.5 * (d1 + d2);
    d2_d1 = 1.0 / d1_d2;

    pdfCoef = B.getInverseBetaFunction();
    pdfCoef *= std::pow(d1_d2, a + 1);
}

void FisherSnedecorRand::setFirstDegree(int degree1)
{
    d1 = std::max(degree1, 1);

    B.setAlpha(.5 * d1);

    a = .5 * d1 - 1;
    d1_d2 = static_cast<double>(d1) / d2;
    c = -.5 * (d1 + d2);
    d2_d1 = 1.0 / d1_d2;

    pdfCoef = B.getInverseBetaFunction();
    pdfCoef *= std::pow(d1_d2, a + 1);
}

void FisherSnedecorRand::setSecondDegree(int degree2)
{
    d2 = std::max(degree2, 1);

    B.setBeta(.5 * d2);

    d1_d2 = static_cast<double>(d1) / d2;
    d2_d1 = 1.0 / d1_d2;
    c = -.5 * (d1 + d2);

    pdfCoef = B.getInverseBetaFunction();
    pdfCoef *= std::pow(d1_d2, a + 1);
}

double FisherSnedecorRand::f(double x) const
{
    return (x < 0) ? 0 : pdfCoef * std::pow(x, a) * std::pow(1 + d1_d2 * x, c);
}

double FisherSnedecorRand::F(double x) const
{
    return B.F(d1_d2 * x);
}

double FisherSnedecorRand::variate() const
{
    return d2_d1 * B.variate();
}

void FisherSnedecorRand::sample(QVector<double> &outputData) const
{
    B.sample(outputData);
    for (double &var : outputData)
        var = d2_d1 * var;
}

double FisherSnedecorRand::Mean() const
{
    return (d2 > 2) ? d2 / (d2 - 2) : INFINITY /*or NAN*/;
}

double FisherSnedecorRand::Variance() const
{
    if (d2 <= 4)
        return INFINITY; /*or NAN*/
    double numerator = 2 * d2 * d2 * (d1 + d2 - 2);
    double denominator = d2 - 2;
    denominator *= denominator;
    denominator *= d1 * (d2 - 4);
    return numerator / denominator;
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
