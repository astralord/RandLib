#ifndef EXPONENTIALNORMALRAND_H
#define EXPONENTIALNORMALRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"
#include "ExponentialRand.h"

class RANDLIBSHARED_EXPORT ExponentialNormalRand : public ContinuousDistribution
{
    NormalRand X;
    ExponentialRand Y;

    double a, b, c, v; /// auxiliary variables

public:
    explicit ExponentialNormalRand(double location = 0, double variance = 1, double rate = 1);
    std::string name() override;

    void setParameters(double location, double variance, double rate);
    inline double getLocation() { return X.getLocation(); }
    inline double getScale() { return X.getScale(); }
    inline double getRate() { return Y.getRate(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double rootVar, double rate);

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Skewness() const override;
    double ExcessKurtosis() const override;
};



#endif // EXPONENTIALNORMALRAND_H
