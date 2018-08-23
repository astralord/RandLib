#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include "BetaRand.h"
#include "ParetoRand.h"

/**
 * @brief The UniformRand class <BR>
 * Uniform continuous distribution
 *
 * f(x | a, b) = 1 / (b - a) for a < x < b
 *
 * Notation: X ~ U(a, b)
 *
 * Related distributions: <BR>
 * X ~ B(1, 1, a, b) <BR>
 * (X - a) / (b - a) ~ IH(1)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT UniformRand : public BetaDistribution<RealType>
{
public:
    UniformRand(double minValue = 0, double maxValue = 1);
    String Name() const override;

    using BetaDistribution<RealType>::SetSupport;

    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return this->a; }
    RealType MaxValue() const override { return this->b; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    /**
     * @fn StandardVariate
     * @param randGenerator
     * @return a random number on interval (0,1) if no preprocessors are specified
     */
    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

    /**
     * @fn StandardVariateClosed
     * @param randGenerator
     * @return a random number on interval [0,1]
     */
    static RealType StandardVariateClosed(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

    /**
     * @fn StandardVariateHalfClosed
     * @param randGenerator
     * @return a random number on interval [0,1)
     */
    static RealType StandardVariateHalfClosed(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

    void Sample(std::vector<RealType> &outputData) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

    static constexpr char TOO_LARGE_A[] = "Minimum element of the sample is smaller than lower boundary returned by method: ";
    static constexpr char TOO_SMALL_B[] = "Maximum element of the sample is greater than upper boundary returned by method: ";

public:
    inline double Entropy() const;

    double LikelihoodFunction(const std::vector<RealType> &sample) const override;
    double LogLikelihoodFunction(const std::vector<RealType> &sample) const override;

    /**
     * @fn FitMinimum
     * fit minimum with maximum-likelihood estimator if unbiased == false,
     * fit minimum using UMVU estimator otherwise
     * @param sample
     */
    void FitMinimum(const std::vector<RealType> &sample, bool unbiased = false);
    /**
     * @fn FitMaximum
     * fit maximum with maximum-likelihood estimator if unbiased == false,
     * fit maximum using UMVU estimator otherwise
     * @param sample
     */
    void FitMaximum(const std::vector<RealType> &sample, bool unbiased = false);
    /**
     * @fn Fit
     * fit support with maximum-likelihood estimator if unbiased == false,
     * fit support using UMVU estimator otherwise
     * @param sample
     */
    void Fit(const std::vector<RealType> &sample, bool unbiased = false);

    /**
     * @fn FitMaximumBayes
     * fit maximum parameter, using Bayesian estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    ParetoRand<RealType> FitMaximumBayes(const std::vector<RealType> &sample, const ParetoRand<RealType> &priorDistribution, bool MAP = false);
};

#endif // UNIFORMRAND_H
