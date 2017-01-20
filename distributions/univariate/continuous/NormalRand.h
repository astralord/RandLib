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
    static double stairWidth[257], stairHeight[256];
    static constexpr double x1 = 3.6541528853610088;
    static const bool dummy;
    static bool SetupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    std::string Name() const override;

private:
    using StableRand::SetParameters;
    using StableRand::GetExponent;
    using StableRand::GetSkewness;

public:
    void SetScale(double scale);
    void SetVariance(double var);
    inline double GetScale() const { return sigma0; }
    inline double GetPrecision() const { return 1.0 / (sigma0 * sigma0); }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    static double Variate(double mean, double rootVar);
    static double StandardVariate();

    static double standardQuantile(double p);
    static double quantile(double p, double mean, double scale);

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Moment(int n) const;

    /// Maximum likelihood estimators
    bool FitMeanMLE(const std::vector<double> &sample);
    bool FitVarianceMLE(const std::vector<double> &sample);
    bool FitMLE(const std::vector<double> &sample);

    /// Method of moments (the same as MLE)
    bool FitMeanMM(const std::vector<double> &sample);
    bool FitVarianceMM(const std::vector<double> &sample);
    bool FitMM(const std::vector<double> &sample);
    
    /// Uniformly minimum variance unbiased (UMVU) estimators
    bool FitMeanUMVU(const std::vector<double> &sample);
    bool FitVarianceUMVU(const std::vector<double> &sample);
    bool FitUMVU(const std::vector<double> &sample);
    bool FitUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha);

    /// Bayesian estimation
    bool FitMeanBayes(const std::vector<double> &sample, NormalRand &priorDistribution);
    bool FitVarianceBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution);
    bool FitBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
