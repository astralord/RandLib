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
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return INT_MAX; }

    void SetRate(double rate);
    inline double GetRate() const { return lambda; }

    double P(int k) const override;
    double F(int k) const override;
    int Variate() const override;
    static int Variate(double rate);
    void Sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;

    bool FitMLE(const std::vector<int> &Sample);
    bool FitMM(const std::vector<int> &Sample);
    bool FitBayes(const std::vector<int> &Sample, GammaRand & priorDistribution);
};

#endif // POISSONRAND_H
