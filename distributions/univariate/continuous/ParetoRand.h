#ifndef PARETORAND_H
#define PARETORAND_H

#include "ExponentialRand.h"

/**
 * @brief The ParetoRand class <BR>
 * Pareto distribution
 *
 * Notation: X ~ Pareto(α, σ)
 *
 * Related distributions: <BR>
 * ln(X / σ) ~ Exp(α)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< shape α
    double sigma = 1; ///< scale σ
    double logAlpha = 0; ///< log(α)
    double logSigma = 0; ///< log(σ)

public:
    ParetoRand(double shape = 1, double scale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return sigma; }
    RealType MaxValue() const override { return INFINITY; }

    void SetShape(double shape);
    void SetScale(double scale);
    inline double GetShape() const { return alpha; }
    inline double GetScale() const { return sigma; }
    inline double GetLogShape() const { return logAlpha; }
    inline double GetLogScale() const { return logSigma; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

private:
    static RealType variateForAlphaEqualOne(RandGenerator &randGenerator);
    static RealType variateForAlphaEqualTwo(RandGenerator &randGenerator);
    static RealType variateForGeneralAlpha(double shape, RandGenerator &randGenerator);

public:
    RealType Variate() const override;
    static RealType StandardVariate(double shape, RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);
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

public:
    long double Entropy() const override;

    /**
     * @fn FitShape
     * Fit α to maximum-likelihood estimator
     * or to UMVU if unbiased == true
     * @param sample
     * @param unbiased
     */
    void FitShape(const std::vector<RealType> &sample, bool unbiased = false);

    /**
     * @fn FitScale
     * Fit σ to maximum-likelihood estimator
     * or to UMVU if unbiased == true
     * @param sample
     * @param unbiased
     */
    void FitScale(const std::vector<RealType> &sample, bool unbiased = false);

    /**
     * @fn Fit
     * Fit parameters to maximum-likelihood estimators
     * or to UMVU if unbiased == true
     * @param sample
     * @param unbiased
     */
    void Fit(const std::vector<RealType> &sample, bool unbiased = false);

    /**
     * @fn FitShapeBayes
     * Fit α, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution of α
     */
    GammaRand<RealType> FitShapeBayes(const std::vector<RealType> &sample, const GammaDistribution<RealType> &priorDistribution, bool MAP = false);
};

#endif // PARETORAND_H
