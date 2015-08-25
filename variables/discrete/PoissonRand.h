#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"

/**
 * @brief The PoissonRand class
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteRand<int>
{
    double lambda;
    double expLambda; /// exp(-l)
public:
    PoissonRand(double rate);
    virtual std::string name() override;

    void setRate(double rate);
    inline double getRate() const { return lambda; }

    double P(int k) const override;
    double F(double x) const override;
    std::complex<double> CF(double t) const override;
    double variate() const override;

    double E() const override { return lambda; }
    double Var() const override { return lambda; }

    inline double Mode() const { return std::floor(lambda); }
    inline double Skewness() const { return 1.0 / std::sqrt(lambda); }
    inline double ExcessiveKurtosis() const { return 1.0 / lambda; }
};

#endif // POISSONRAND_H
