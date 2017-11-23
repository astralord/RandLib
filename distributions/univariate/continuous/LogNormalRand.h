#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"

/**
 * @brief The LogNormalRand class <BR>
 * Log-Normal distribution
 *
 * Notation X ~ Log-Normal(μ, σ)
 *
 * Related distributions: <BR>
 * ln(X) ~ Normal(μ, σ)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousDistribution<RealType>
{
    NormalRand<RealType> X{};
    double expMu = 1; ///< exp(μ)
    double expHalfSigmaSq = 1.6487212707; ///< exp(σ^2 / 2)

public:
    LogNormalRand(double location = 0, double squaredScale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return X.Mean(); }
    inline double GetScale() const { return X.GetScale(); }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

    RealType Variate() const override;
    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    void FitLocation(const std::vector<RealType> &sample);
    void FitScale(const std::vector<RealType> &sample);
    void Fit(const std::vector<RealType> &sample);

    /**
     * @fn FitLocationBayes
     * Set location, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    NormalRand<RealType> FitLocationBayes(const std::vector<RealType> &sample, const NormalRand<RealType> &priorDistribution, bool MAP = false);

    /**
     * @fn FitScaleBayes
     * Set scale, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    InverseGammaRand<RealType> FitScaleBayes(const std::vector<RealType> &sample, const InverseGammaRand<RealType> &priorDistribution, bool MAP = false);

    /**
     * @fn FitBayes
     * Set parameters, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    NormalInverseGammaRand<RealType> FitBayes(const std::vector<RealType> &sample, const NormalInverseGammaRand<RealType> &priorDistribution, bool MAP = false);
};

#endif // LOGNORMALRAND_H
