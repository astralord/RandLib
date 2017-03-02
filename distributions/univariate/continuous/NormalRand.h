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
 */
class RANDLIBSHARED_EXPORT NormalRand : public StableRand
{
    double sigma0;

    /// Tables for ziggurat
    static long double stairWidth[257], stairHeight[256];
    static constexpr long double x1 = 3.6541528853610088l;
    static const bool dummy;
    static bool SetupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    std::string Name() const override;

private:
    using StableRand::SetParameters;

public:
    void SetScale(double scale);
    void SetVariance(double var);
    inline double GetScale() const { return sigma0; }
    inline double GetPrecision() const { return 1.0 / (sigma0 * sigma0); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    static double StandardVariate();

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Moment(int n) const;

    /// Maximum likelihood estimators
    void FitMeanMLE(const std::vector<double> &sample);
    void FitVarianceMLE(const std::vector<double> &sample);
    void FitMeanAndVarianceMLE(const std::vector<double> &sample);
    
    /// Uniformly minimum variance unbiased (UMVU) estimators
    void FitMeanUMVU(const std::vector<double> &sample);
    void FitVarianceUMVU(const std::vector<double> &sample);
    void FitMeanAndVarianceUMVU(const std::vector<double> &sample);
    void FitMeanAndVarianceUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha);

    /// Bayesian estimation
    NormalRand FitMeanBayes(const std::vector<double> &sample, const NormalRand &priorDistribution);
    InverseGammaRand FitVarianceBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution);
    NormalInverseGammaRand FitMeanAndVarianceBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
