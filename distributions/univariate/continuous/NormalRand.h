#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "StableRand.h"
#include "InverseGammaRand.h"
#include "../../bivariate/NormalInverseGammaRand.h"

/**
 * @brief The NormalRand class
 * Normal distribution
 *
 * f(x | μ, σ) = 1 / ((2 π σ^2)^(1/2) * exp(-(x - μ)^2 / (2 σ^2))
 *
 * Notation: X ~ N(μ, σ)
 *
 * Related distributions:
 * X ~ S(2, 0, √2σ, μ)
 */
class RANDLIBSHARED_EXPORT NormalRand : public StableDistribution
{
    double sigma;

    /// Tables for ziggurat
    static long double stairWidth[257], stairHeight[256];
    static constexpr long double x1 = 3.6541528853610088l;
    static const bool dummy;
    static bool SetupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    std::string Name() const override;

public:
    void SetScale(double scale);
    void SetVariance(double var);
    /**
     * @brief GetScale
     * @return σ
     */
    inline double GetScale() const { return sigma; }
    /**
     * @brief GetLogScale
     * @return log(σ)
     */
    inline double GetLogScale() const { return logGamma - 0.5 * M_LN2; }
    /**
     * @brief GetPrecision
     * @return 1/σ^2
     */
    inline double GetPrecision() const { return 1.0 / (sigma * sigma); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    /**
     * @brief StandardVariate
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

    /**
     * @brief FitMeanMLE
     * set mean, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitMeanMLE(const std::vector<double> &sample);
    /**
     * @brief FitVarianceMLE
     * set variance, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitVarianceMLE(const std::vector<double> &sample);
    /**
     * @brief FitMeanAndVarianceMLE
     * set mean and variance, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitMeanAndVarianceMLE(const std::vector<double> &sample);
    /**
     * @brief FitMeanUMVU
     * set mean, returned by uniformly minimum variance unbiased estimator
     * @param sample
     */
    void FitMeanUMVU(const std::vector<double> &sample);
    /**
     * @brief FitVarianceUMVU
     * set variance, returned by uniformly minimum variance unbiased estimator
     * @param sample
     */
    void FitVarianceUMVU(const std::vector<double> &sample);
    /**
     * @brief FitMeanAndVarianceUMVU
     * set mean and variance, returned by uniformly minimum variance unbiased estimator
     * @param sample
     */
    void FitMeanAndVarianceUMVU(const std::vector<double> &sample);
    /**
     * @brief FitMeanAndVarianceUMVU
     * set mean and variance, returned by uniformly minimum variance unbiased estimator
     * and returns confidence intervals for these parameters
     * @param sample
     * @param confidenceIntervalForMean
     * @param confidenceIntervalForVariance
     * @param alpha
     */
    void FitMeanAndVarianceUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha);
    /**
     * @brief FitMeanBayes
     * set mean, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    NormalRand FitMeanBayes(const std::vector<double> &sample, const NormalRand &priorDistribution);
    /**
     * @brief FitVarianceBayes
     * set variance, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    InverseGammaRand FitVarianceBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution);
    /**
     * @brief FitMeanAndVarianceBayes
     * set mean and variance, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    NormalInverseGammaRand FitMeanAndVarianceBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
