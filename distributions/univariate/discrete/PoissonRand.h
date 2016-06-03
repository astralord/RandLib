#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/GammaRand.h"

/**
 * @brief The PoissonRand class
 *
 * P(X = k) = \lambda^k * e^(-\lambda) / k!
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteDistribution
{
    double lambda;
    double expmLambda; /// exp(-lambda)
    double logLambda; /// ln(lambda)
    int floorLambda; /// floor(lambda)
    double FFloorLambda; /// P(X < floor(lambda))
    double PFloorLambda; /// P(X = floor(lambda))

public:
    explicit PoissonRand(double rate = 1.0);
    std::string name() override;

    void setRate(double rate);
    inline double getRate() const { return lambda; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;
    static int variate(double rate);

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool checkValidity(const std::vector<double> &sample);

    bool fitRateMLE(const std::vector<double> &sample);
    bool fitRateBayes(const std::vector<double> &sample, GammaRand & priorDistribution);
};

#endif // POISSONRAND_H
