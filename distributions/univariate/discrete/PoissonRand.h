#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/GammaRand.h"

/**
 * @brief The PoissonRand class
 * Poisson distribution
 *
 * P(X = k) = λ^k * e^(-λ) / k!
 *
 * Notation: X ~ Po(λ)
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteDistribution
{
    double lambda;
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

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    static int Variate(double rate);

    double Mean() const override;
    double Variance() const override;
    int Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    void FitRateMLE(const std::vector<int> &sample);
    void FitRateMM(const std::vector<int> &sample);
    GammaRand FitRateBayes(const std::vector<int> &sample, const GammaDistribution & priorDistribution);
};

#endif // POISSONRAND_H
