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
class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousDistribution
{
    NormalRand X{};
    double expMu = 1; ///< exp(μ)
    double expHalfSigmaSq = 1.6487212707; ///< exp(σ^2 / 2)

public:
    LogNormalRand(double location = 0, double squaredScale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return X.Mean(); }
    inline double GetScale() const { return X.GetScale(); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;

    double Variate() const override;
    static double StandardVariate(RandGenerator &randGenerator = staticRandGenerator);
    void Reseed(unsigned long seed) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    void FitLocation(const std::vector<double> &sample);
    void FitScale(const std::vector<double> &sample);
    void Fit(const std::vector<double> &sample);

    /**
     * @fn FitLocationBayes
     * Set location, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    NormalRand FitLocationBayes(const std::vector<double> &sample, const NormalRand &priorDistribution, bool MAP = false);

    /**
     * @fn FitScaleBayes
     * Set scale, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    InverseGammaRand FitScaleBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution, bool MAP = false);

    /**
     * @fn FitBayes
     * Set parameters, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    NormalInverseGammaRand FitBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution, bool MAP = false);
};

#endif // LOGNORMALRAND_H
