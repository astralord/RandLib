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
template < typename IntType = int >
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteDistribution<IntType>
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
    IntType MinValue() const override { return 0; }
    IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

private:
    void SetGeneratorConstants();

public:
    void SetRate(double rate);
    inline double GetRate() const { return lambda; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;
private:
    double acceptanceFunction(IntType X) const;
    bool generateByInversion() const;
    IntType variateRejection() const;
    IntType variateInversion() const;

public:
    IntType Variate() const override;
    static IntType Variate(double rate, RandGenerator &randGenerator = ProbabilityDistribution<IntType>::staticRandGenerator);
    void Sample(std::vector<IntType> &outputData) const;

    long double Mean() const override;
    long double Variance() const override;
    IntType Median() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @fn Fit
     * fit rate λ via maximum-likelihood method
     * @param sample
     */
    void Fit(const std::vector<IntType> &sample);
    /**
     * @brief Fit
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     */
    void Fit(const std::vector<IntType> &sample, DoublePair &confidenceInterval, double significanceLevel);
    /**
     * @fn FitBayes
     * fit rate λ via Bayes estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior Gamma distribution
     */
    GammaRand<> FitBayes(const std::vector<IntType> &sample, const GammaDistribution<> & priorDistribution, bool MAP = false);
};

#endif // POISSONRAND_H
