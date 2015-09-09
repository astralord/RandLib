#ifndef BERNOULLIRAND_H
#define BERNOULLIRAND_H

#include "DiscreteRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The BernoulliRand class
 */
class RANDLIBSHARED_EXPORT BernoulliRand : public DiscreteRand<int>
{
    double p, q;
    unsigned long long generatorEdge; /// such value that probability of (BasicRandGenerator's variate > generatorEdge) is equal p

public:
    explicit BernoulliRand(double probability = 0.5);
    virtual std::string name() override;

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double p);

    double E() const override { return p; }
    double Var() const override { return p * (1 - p); }

    std::complex<double> CF(double t) const override;

    inline double Median();
    inline double Skewness();
    inline double ExcessiveKurtosis();

    inline double Entropy();

};

#endif // BERNOULLIRAND_H
