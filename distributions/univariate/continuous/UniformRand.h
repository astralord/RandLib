#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include "BetaRand.h"

/**
 * @brief The UniformRand class
 * Uniform dontinuous distribution
 * X ~ U(a, b)
 *
 * f(x|a, b) = 1 / (b - a) for a < x < b
 *
 * X ~ Beta(1, 1)
 */
class RANDLIBSHARED_EXPORT UniformRand : public BetaRand
{
    double bmaInv; /// 1 / (b - a)

public:
    UniformRand(double minValue = 0, double maxValue = 1);
    std::string name() override;

    void setSupport(double minValue, double maxValue);
    inline double getMinValue() const { return a; }
    inline double getMaxValue() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double minValue, double maxValue);
    static double standardVariate();

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() const;

    /// Maximum likelihood estimation
    bool fitMinimumMLE(const QVector<double> &sample);
    bool fitMaximumMLE(const QVector<double> &sample);
    bool fitSupportMLE(const QVector<double> &sample);
    
    /// Method of moments
    bool fitMinimumMM(const QVector<double> &sample);
    bool fitMaximumMM(const QVector<double> &sample);
    bool fitSupportMM(const QVector<double> &sample);
    
    /// Minimum-variance unbiased estimator
    bool fitMinimumUMVU(const QVector<double> &sample);
    bool fitMaximumUMVU(const QVector<double> &sample);
    bool fitSupportUMVU(const QVector<double> &sample);
};

#endif // UNIFORMRAND_H
