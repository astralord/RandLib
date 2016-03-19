#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

/**
 * @brief The PoissonRand class
 *
 * P(X = k) = l^k * e^(-l) / k!
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteDistribution
{
    double lambda;
    double expmLambda; /// exp(-lambda)
    double logLambda; /// ln(lambda)
    int floorLambda; /// floor(lambda) and
    double FFloorLambda; /// P(X < floor(lambda))
    double PFloorLambda; /// P(X = floor(lambda))

public:
    explicit PoissonRand(double rate = 1.0);
    std::string name() override;

    void setRate(double rate);
    inline double getRate() const { return lambda; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double rate);

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool fit_MLE(const QVector<int> &sample);
};

#endif // POISSONRAND_H
