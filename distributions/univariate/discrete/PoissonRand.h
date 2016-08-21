#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/GammaRand.h"

/**
 * @brief The PoissonRand class
 *
 * P(X = k) = λ^k * e^(-λ) / k!
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteDistribution
{
    double lambda;
    double expmLambda; /// exp(-λ)
    double logLambda; /// ln(λ)
    int floorLambda; /// [λ]
    double FFloorLambda; /// P(X < [λ])
    double PFloorLambda; /// P(X = [λ])

public:
    explicit PoissonRand(double rate = 1.0);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return INT_MAX; }

    void setRate(double rate);
    inline double getRate() const { return lambda; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;
    static int variate(double rate);
    void sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;

    bool fitMLE(const std::vector<int> &sample);
    bool fitMM(const std::vector<int> &sample);
    bool fitBayes(const std::vector<int> &sample, GammaRand & priorDistribution);
};

#endif // POISSONRAND_H
