#ifndef CHISQUAREDRAND_H
#define CHISQUAREDRAND_H

#include <RandomVariable.h>
#include <NormalRand.h>

class RANDLIBSHARED_EXPORT ChiSquaredRand : public RandomVariable
{
    int k;
    double pdfCoef;
    NormalRand X;

    // TODO: create math class for it
    double factorial(int n);
    double doubleFactorial(int n);

public:
    ChiSquaredRand(int k);

    void setDegree(int degrees);
    size_t getDegree() { return k; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return (double)k; }
    inline double Var() { return k + k; }

    inline double Mode() { return std::max(k - 2, 0); }
    inline double Skewness() { return std::sqrt(8.0 / k); }
    inline double ExcessKurtosis() { return 12.0 / k; }
};

#endif // CHISQUAREDRAND_H
