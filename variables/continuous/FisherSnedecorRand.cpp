#include "FisherSnedecorRand.h"

FisherSnedecorRand::FisherSnedecorRand(int degree1, int degree2)
{
    setDegrees(degree1, degree2);
}

void FisherSnedecorRand::setDegrees(int degree1, int degree2)
{
    d1 = std::max(degree1, 1);
    U1.setDegree(d1);
    gammaA = RandMath::gammaHalf(d1);

    d2 = std::max(degree2, 1);
    U2.setDegree(d2);
    gammaB = RandMath::gammaHalf(d2);

    a = .5 * d1 - 1;
    b = static_cast<double>(d1) / d2;
    c = -.5 * (d1 + d2);

    pdfCoef = RandMath::gammaHalf(d1 + d2) / (gammaA * gammaB);
    pdfCoef *= std::pow(b, a + 1);
}

void FisherSnedecorRand::setFirstDegree(int degree1)
{
    d1 = std::max(degree1, 1);
    U1.setDegree(d1);
    gammaA = RandMath::gammaHalf(d1);

    a = .5 * d1 - 1;
    b = static_cast<double>(d1) / d2;
    c = -.5 * (d1 + d2);

    pdfCoef = RandMath::gammaHalf(d1 + d2) / (gammaA * gammaB);
    pdfCoef *= std::pow(b, a + 1);
}

void FisherSnedecorRand::setSecondDegree(int degree2)
{
    d2 = std::max(degree2, 1);
    U2.setDegree(d2);
    gammaB = RandMath::gammaHalf(d2);

    b = static_cast<double>(d1) / d2;
    c = -.5 * (d1 + d2);

    pdfCoef = RandMath::gammaHalf(d1 + d2) / (gammaA * gammaB);
    pdfCoef *= std::pow(b, a + 1);
}

double FisherSnedecorRand::f(double x) const
{
    return (x >= 0) ? pdfCoef * std::pow(x, a) * std::pow(1 + b * x, c) : 0;
}

double FisherSnedecorRand::F(double x) const
{
    return x;
}

double FisherSnedecorRand::value()
{
    double numerator = d2 * U1.value();
    double denominator = d1 * U2.value();
    return numerator / denominator;
}
