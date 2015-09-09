#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

/**
 * @brief The PoissonRand class
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteRand<int>
{
    double lambda;
    double expLambda; /// exp(-l)
public:
    explicit PoissonRand(double rate);
    virtual std::string name() override;

    void setRate(double rate);
    inline double getRate() const { return lambda; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double rate);

    double E() const override { return lambda; }
    double Var() const override { return lambda; }

    std::complex<double> CF(double t) const override;

    inline double Mode() const;
    inline double Skewness() const;
    inline double ExcessiveKurtosis() const;
};

#endif // POISSONRAND_H
