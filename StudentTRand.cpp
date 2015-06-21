#include "StudentTRand.h"

StudentTRand::StudentTRand(int degree) :
    X(0, 1),
    Y(degree)
{
    setDegree(degree);
}

void StudentTRand::setDegree(int degree)
{
    v = degree;
    Y.setDegree(v);
}

double StudentTRand::pdf(double x)
{
    return x;
}

double StudentTRand::cdf(double x)
{
    return x;
}

double StudentTRand::value()
{
    //v = 1 -> cauchy
    //v = inf -> normal
    return X.value() * qSqrt((double)v / Y.value());
}
