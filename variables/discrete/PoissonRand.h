#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

/**
 * @brief The PoissonRand class
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteRand
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

    double Mean() const override { return lambda; }
    double Variance() const override { return lambda; }

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // POISSONRAND_H
