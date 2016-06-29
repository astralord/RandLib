#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "StableRand.h"
#include "InverseGammaRand.h"
#include "../../multivariate/NormalInverseGammaRand.h"

/**
 * @brief The NormalRand class
 * Normal distribution
 * X ~ N(mu, sigma)
 *
 * f(x|\mu, \sigma) = 1 / (\sqrt(2 \pi) \sigma) * exp(-(x - \mu)^2 / (2 \sigma^2))
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
    double Quantile(double p) const override;

    double Moment(int n) const;

    /// Maximum likelihood estimators
    bool fitMeanMLE(const std::vector<double> &sample);
    bool fitVarianceMLE(const std::vector<double> &sample);
    bool fitMeanAndVarianceMLE(const std::vector<double> &sample);

    /// Method of moments (actually, the same as MLE)
    bool fitMeanMM(const std::vector<double> &sample);
    bool fitVarianceMM(const std::vector<double> &sample);
    bool fitMeanAndVarianceMM(const std::vector<double> &sample);
    
    /// Uniformly minimum variance unbiased (UMVU) estimators
    bool fitMeanUMVU(const std::vector<double> &sample);
    bool fitVarianceUMVU(const std::vector<double> &sample);
    bool fitMeanAndVarianceUMVU(const std::vector<double> &sample);
    bool fitMeanAndVarianceUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha);

    /// Bayesian estimation
    bool fitMeanBayes(const std::vector<double> &sample, NormalRand &priorDistribution);
    bool fitVarianceBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution);
    bool fitMeanAndVarianceBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
