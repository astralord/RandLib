#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/GammaRand.h"

/**
 * @brief The PoissonRand class <BR>
 * Poisson distribution
 *
 * P(X = k) = λ^k * exp(-λ) / k!
 *
 * Notation: X ~ Po(λ)
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteDistribution
{
    double lambda = 1; ///< rate λ
    double logLambda = 0; ///< ln(λ)
    double Fmu = 2 * M_1_E; ///< P(X < [λ])
    double Pmu = M_1_E; ///< P(X = [λ])

    double mu = 1, delta = 6;
    double zeta{};
    long double c1{}, c2{}, c3{}, c4{}, c{};
    double sqrtMu = 1, sqrtMupHalfDelta = 2;
    double lfactMu = 0;

public:
    explicit PoissonRand(double rate = 1.0);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return INT_MAX; }

private:
    void SetGeneratorConstants();

public:
    void SetRate(double rate);
    inline double GetRate() const { return lambda; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
private:
    double acceptanceFunction(int X) const;
    bool generateByInversion() const;
    int variateRejection() const;
    int variateInversion() const;

public:
    int Variate() const override;
    static int Variate(double rate);
    void Sample(std::vector<int> &outputData) const;

    double Mean() const override;
    double Variance() const override;
    int Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @fn Fit
     * fit rate λ via maximum-likelihood method
     * @param sample
     */
    void Fit(const std::vector<int> &sample);
    /**
     * @brief Fit
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     */
    void Fit(const std::vector<int> &sample, DoublePair &confidenceInterval, double significanceLevel);
    /**
     * @fn FitBayes
     * fit rate λ via Bayes estimation
     * @param sample
     * @param priorDistribution
     * @return posterior Gamma distribution
     */
    GammaRand FitBayes(const std::vector<int> &sample, const GammaDistribution & priorDistribution);
};

#endif // POISSONRAND_H
