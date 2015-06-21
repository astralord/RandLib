#include "FisherSnedecorRand.h"

FisherSnedecorRand::FisherSnedecorRand(int degree1, int degree2) :
    U1(degree1),
    U2(degree2)
{
    setDegrees(degree1, degree2);
}

void FisherSnedecorRand::setDegrees(int degree1, int degree2)
{
    d1 = degree1;
    d2 = degree2;
    U1.setDegree(d1);
    U2.setDegree(d2);
}

double FisherSnedecorRand::pdf(double x)
{
    return x;
}

double FisherSnedecorRand::cdf(double x)
{
    return x;
}

double FisherSnedecorRand::value()
{
    double numen = d2 * U1.value();
    double denom = d1 * U2.value();
    return numen / denom;
}
