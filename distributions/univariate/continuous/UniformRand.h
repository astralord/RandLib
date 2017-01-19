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

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

private:
    using BetaRand::SetParameters;

public:
    void SetSupport(double minValue, double maxValue);
    inline double GetMinValue() const { return a; }
    inline double GetMaxValue() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    static double Variate(double minValue, double maxValue);
    static double StandardVariate();

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    inline double Entropy() const;

    /// Maximum likelihood estimation
    bool FitMinimumMLE(const std::vector<double> &sample);
    bool FitMaximumMLE(const std::vector<double> &sample);
    bool FitMLE(const std::vector<double> &sample);
    
    /// Method of moments
    bool FitMinimumMM(const std::vector<double> &sample);
    bool FitMaximumMM(const std::vector<double> &sample);
    bool FitMM(const std::vector<double> &sample);
    
    /// Minimum-variance unbiased estimator
    bool FitMinimumUMVU(const std::vector<double> &sample);
    bool FitMaximumUMVU(const std::vector<double> &sample);
    bool FitUMVU(const std::vector<double> &sample);
};

#endif // UNIFORMRAND_H
