#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "StableRand.h"
#include "InverseGammaRand.h"
#include "../../bivariate/NormalInverseGammaRand.h"

/**
 * @brief The NormalRand class <BR>
 * Normal distribution
 *
 * f(x | μ, σ^2) = 1 / ((2 π σ^2)^(1/2) * exp(-(x - μ)^2 / (2 σ^2))
 *
 * Notation: X ~ N(μ, σ^2)
 *
 * Related distributions: <BR>
 * X ~ S(2, 0, σ/√2, μ)
 */
class RANDLIBSHARED_EXPORT NormalRand : public StableDistribution
{
    double sigma = 1; ///< scale σ

    static long double stairWidth[257]; ///< width of ziggurat's stairs
    static long double stairHeight[256]; ///< height of ziggurat's stairs
    static constexpr long double x1 = 3.6541528853610088l; ///< starting point for ziggurat's setup
    static const bool dummy;
    static bool SetupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    String Name() const override;

public:
    void SetScale(double scale);
    void SetVariance(double var);
    /**
     * @fn GetScale
     * @return σ
     */
    inline double GetScale() const { return sigma; }
    /**
     * @fn GetLogScale
     * @return log(σ)
     */
    inline double GetLogScale() const { return StableDistribution::GetLogScale() - 0.5 * M_LN2; }
    /**
     * @fn GetPrecision
     * @return 1/σ^2
     */
    inline double GetPrecision() const { return 1.0 / (sigma * sigma); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    /**
     * @fn StandardVariate
     * @return variate from standard normal distribution
     */
    static double StandardVariate();
    void Sample(std::vector<double> &outputData) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Moment(int n) const;
    double ThirdMoment() const override { return Moment(3); }
    double FourthMoment() const override { return Moment(4); }

    /**
     * @fn FitLocation
     * set location, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitLocation(const std::vector<double> &sample);

    /**
     * @brief FitLocation
     * set location, returned by maximium-likelihood estimator
     * and return confidenceInterval for given significance level
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     */
    void FitLocation(const std::vector<double> &sample, DoublePair &confidenceInterval, double significanceLevel);

    /**
     * @fn FitVariance
     * set variance, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitVariance(const std::vector<double> &sample);

    /**
     * @brief FitVariance
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     * @param unbiased
     */
    void FitVariance(const std::vector<double> &sample, DoublePair &confidenceInterval, double significanceLevel, bool unbiased = false);

    /**
     * @brief FitScale
     * set scale, returned via maximum-likelihood estimation or unbiased estimator
     * (which might be different from the unbiased estimator of variance)
     * @param sample
     * @param unbiased
     */
    void FitScale(const std::vector<double> &sample, bool unbiased = false);

    /**
     * @fn Fit
     * set parameters, returned by maximium-likelihood estimator if unbiased = false,
     * otherwise set parameters via UMVU estimator
     * @param sample
     * @param unbiased
     */
    void Fit(const std::vector<double> &sample, bool unbiased = false);

    /**
     * @brief Fit
     * set parameters, returned by maximium-likelihood estimator if unbiased = false,
     * otherwise set parameters via UMVU estimator, and return confidence intervals for given significance level
     * @param sample
     * @param confidenceIntervalForMean
     * @param confidenceIntervalForVariance
     * @param significanceLevel
     * @param unbiased
     */
    void Fit(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double significanceLevel, bool unbiased = false);

    /**
     * @fn FitLocationBayes
     * set location, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    NormalRand FitLocationBayes(const std::vector<double> &sample, const NormalRand &priorDistribution);

    /**
     * @fn FitVarianceBayes
     * set variance, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    InverseGammaRand FitVarianceBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution);

    /**
     * @fn FitBayes
     * set parameters, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    NormalInverseGammaRand FitBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
