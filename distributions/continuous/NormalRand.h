#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "StableRand.h"
#include "InverseGammaRand.h"

/**
 * @brief The NormalRand class
 * Normal distribution
 * X ~ N(mu, sigma)
 *
 * f(x|\mu, \sigma) = 1 / (\sqrt(2 \pi) \sigma) * exp(-(x - \mu)^2 / (2 \sigma^2))
 */
class RANDLIBSHARED_EXPORT NormalRand : public StableRand
{
    //TODO: find a way to initialize them without dummy
    /// Tables for ziggurat
    static double stairWidth[257], stairHeight[256];
    static constexpr double x1 = 3.6541528853610088;
    static const bool dummy;
    static bool setupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    std::string name() override;

private:
    using StableRand::setParameters;

public:
    void setVariance(double var);
    inline double getVariance() const { return sigma * sigma; }
    inline double getPrecision() const { return 1.0 / (sigma * sigma); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double mean, double rootVar);
    static double standardVariate();

    double Mean() const;
    double Variance() const;

    std::complex<double> CF(double t) const override;
    static double standardQuantile(double p);
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Moment(int n) const;

    /// Maximum likelihood estimators
    bool fitMeanMLE(const QVector<double> &sample);
    bool fitVarianceMLE(const QVector<double> &sample);
    bool fitMeanAndVarianceMLE(const QVector<double> &sample);

    /// Method of moments
    bool fitMeanMM(const QVector<double> &sample);
    bool fitVarianceMM(const QVector<double> &sample);
    bool fitMeanAndVarianceMM(const QVector<double> &sample);
    
    /// Uniformly minimum variance unbiased (UMVU) estimators
    bool fitMeanUMVU(const QVector<double> &sample);
    bool fitVarianceUMVU(const QVector<double> &sample);
    bool fitMeanAndVarianceUMVU(const QVector<double> &sample);

    /// Bayesian estimation
    bool fitMeanBayes(const QVector<double> &sample, NormalRand &priorDistribution);
    bool fitVarianceBayes(const QVector<double> &sample, InverseGammaRand &priorDistribution);
};

#endif // NORMALRAND_H
