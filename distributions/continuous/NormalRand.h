#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The NormalRand class
 *
 * f(x|mu, sigma) = 1 / (\sqrt(2 pi) sigma) * exp(-(x - mu)^2 / (2 sigma^2))
 *
 * Normal distribution: X ~ N(mu, sigma)
 */
class RANDLIBSHARED_EXPORT NormalRand : public ContinuousDistribution
{
    double mu, sigma;
    double sigmaSqrt2Inv; /// 1 / (sigma * sqrt(2))

    //TODO: find a way to initialize them without dummy
    /// Tables for ziggurat
    static double stairWidth[257], stairHeight[256];
    static constexpr double x1 = 3.6541528853610088;
    static const bool dummy;
    /// Set up ziggurat tables
    static bool setupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    std::string name() override;

    void setMean(double mean);
    void setSigma(double rootVar);
    void setVariance(double var);
    inline double getMean() const { return mu; }
    inline double getSigma() const { return sigma; }
    inline double getVar() const { return sigma * sigma; }

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
    bool fitMean_MLE(const QVector<double> &sample);
    bool fitVariance_MLE(const QVector<double> &sample);
    bool fit_MLE(const QVector<double> &sample);
    
    /// Method of moments
    bool fitMean_MM(const QVector<double> &sample);
    bool fitVariance_MM(const QVector<double> &sample);
    bool fit_MM(const QVector<double> &sample);
    
    /// Uniformly minimum variance unbiased (UMVU) estimators
    bool fitMean_UMVU(const QVector<double> &sample);
    bool fitVariance_UMVU(const QVector<double> &sample);
    bool fit_UMVU(const QVector<double> &sample);
};

#endif // NORMALRAND_H
