#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "StableRand.h"
#include "InverseGammaRand.h"
#include "../../multivariate/NormalInverseGammaRand.h"

/**
 * @brief The NormalRand class
 * Normal distribution
 *
 * f(x | μ, σ) = 1 / ((2 π σ^2)^(1/2) * exp(-(x - μ)^2 / (2 σ^2))
 *
 * Notation: X ~ N(μ, σ)
 *
 */
class RANDLIBSHARED_EXPORT NormalRand : public StableRand
{
    double sigma0;

    //TODO: find a way to initialize them without dummy
    // try to make tables also constexpr
    /// Tables for ziggurat
    static double stairWidth[257], stairHeight[256];
    static constexpr double x1 = 3.6541528853610088;
    static const bool dummy;
    static bool setupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    std::string name() const override;

private:
    using StableRand::setParameters;
    using StableRand::getExponent;
    using StableRand::getSkewness;

public:
    void setScale(double scale);
    void setVariance(double var);
    inline double getScale() const { return sigma0; }
    inline double getPrecision() const { return 1.0 / (sigma0 * sigma0); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double mean, double rootVar);
    static double standardVariate();

    std::complex<double> CF(double t) const override;
    static double standardQuantile(double p);
    static double quantile(double p, double mean, double scale);

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double Moment(int n) const;

    /// Maximum likelihood estimators
    bool fitMeanMLE(const std::vector<double> &sample);
    bool fitVarianceMLE(const std::vector<double> &sample);
    bool fitMLE(const std::vector<double> &sample);

    /// Method of moments (the same as MLE)
    bool fitMeanMM(const std::vector<double> &sample);
    bool fitVarianceMM(const std::vector<double> &sample);
    bool fitMM(const std::vector<double> &sample);
    
    /// Uniformly minimum variance unbiased (UMVU) estimators
    bool fitMeanUMVU(const std::vector<double> &sample);
    bool fitVarianceUMVU(const std::vector<double> &sample);
    bool fitUMVU(const std::vector<double> &sample);
    bool fitUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha);

    /// Bayesian estimation
    bool fitMeanBayes(const std::vector<double> &sample, NormalRand &priorDistribution);
    bool fitVarianceBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution);
    bool fitBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
