#ifndef BERNOULLIRAND_H
#define BERNOULLIRAND_H

#include "DiscreteRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The BernoulliRand class
 */
class RANDLIBSHARED_EXPORT BernoulliRand : public DiscreteRand<int>
{
    double p;
    unsigned long generatorEdge; /// such value that probability of (BasicRandGenerator's variate > generatorEdge) is equal p

public:
    BernoulliRand(double successProbability);
    void setSuccessProbability(double successProbability);
    inline double getSuccessProbability() { return p; }

    virtual double P(int k) const override;
    virtual double F(double x) const override;
    virtual int variate() override;

    double E() const override { return p; }
    double Var() const override { return p * (1 - p); }

    inline double Median() { return (p < 0.5) ? 0 : ((p > 0.5) ? 1 : 0.5); }
    inline double Skewness() { return (1 - 2 * p) / std::sqrt(p * (1 - p)); }
    inline double ExcessiveKurtosis() { return 1.0 / (p * (1 - p)) - 6; }

    inline double Entropy() { return -(p * std::log(p) + (1- p) * std::log(1 - p)); }
};

#endif // BERNOULLIRAND_H
