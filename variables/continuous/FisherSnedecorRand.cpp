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

    B.setParameters(.5 * d1, .5 * d2);

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
    if (x <= 0)
        return 0;

    return RandMath::integral([this] (double t)
    {
        return f(t);
    },
    0, x);
}

double FisherSnedecorRand::variate() const
{
    return d2_d1 * B.variate();
}

void FisherSnedecorRand::sample(QVector<double> &outputData)
{
    B.sample(outputData);
    for (double &var : outputData)
        var = d2_d1 * var;
}
