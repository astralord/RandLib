#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include "BetaRand.h"

/**
 * @brief The UniformRand class
 * Uniform continuous distribution
 *
 * f(x | a, b) = 1 / (b - a) for a < x < b
 *
 * Notation: X ~ U(a, b)
 *
 * Related distributions:
 * X ~ Beta(1, 1, a, b)
 * (X - a) / (b - a) ~ IH(1)
 */
class RANDLIBSHARED_EXPORT UniformRand : public BetaRand
{
    double bmaInv; /// 1 / (b - a)

public:
    UniformRand(double minValue = 0, double maxValue = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

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
    bool fitMinimumMLE(const std::vector<double> &sample);
    bool fitMaximumMLE(const std::vector<double> &sample);
    bool fitSupportMLE(const std::vector<double> &sample);
    
    /// Method of moments
    bool fitMinimumMM(const std::vector<double> &sample);
    bool fitMaximumMM(const std::vector<double> &sample);
    bool fitSupportMM(const std::vector<double> &sample);
    
    /// Minimum-variance unbiased estimator
    bool fitMinimumUMVU(const std::vector<double> &sample);
    bool fitMaximumUMVU(const std::vector<double> &sample);
    bool fitSupportUMVU(const std::vector<double> &sample);
};

#endif // UNIFORMRAND_H
