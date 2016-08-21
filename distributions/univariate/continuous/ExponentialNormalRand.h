#ifndef EXPONENTIALNORMALRAND_H
#define EXPONENTIALNORMALRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"
#include "ExponentialRand.h"

/**
 * @brief The ExponentialNormalRand class
 * Exponentially modified Gaussian distribution
 *
 * Notation: X ~ EMG(μ, σ, β)
 *
 * Related distributions:
 * X = Y + Z, where Y ~ Normal(μ, σ) and Z ~ Exp(β)
 */
class RANDLIBSHARED_EXPORT ExponentialNormalRand : public ContinuousDistribution
{
    NormalRand X;
    ExponentialRand Y;

    double a, b, c, v; /// auxiliary variables

public:
    explicit ExponentialNormalRand(double location = 0, double variance = 1, double rate = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void setParameters(double location, double variance, double rate);
    inline double getLocation() const { return X.getLocation(); }
    inline double getScale() const { return X.getScale(); }
    inline double getRate() const { return Y.getRate(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double rootVar, double rate);

    double Mean() const override;
    double Variance() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;
};



#endif // EXPONENTIALNORMALRAND_H
