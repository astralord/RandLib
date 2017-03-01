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
public:
    UniformRand(double minValue = 0, double maxValue = 1);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

private:
    using BetaRand::SetParameters;

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
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

    static constexpr char TOO_LARGE_A[] = "Minimum element of the sample is larger than lower boundary returned by method: ";
    static constexpr char TOO_SMALL_B[] = "Maximum element of the sample is smaller than upper boundary returned by method: ";
public:
    inline double Entropy() const;

    /// Maximum likelihood estimation
    void FitMinimumMLE(const std::vector<double> &sample);
    void FitMaximumMLE(const std::vector<double> &sample);
    void FitMLE(const std::vector<double> &sample);
    
    /// Method of moments
    void FitMinimumMM(const std::vector<double> &sample);
    void FitMaximumMM(const std::vector<double> &sample);
    void FitMM(const std::vector<double> &sample);
    
    /// Minimum-variance unbiased estimator
    void FitMinimumUMVU(const std::vector<double> &sample);
    void FitMaximumUMVU(const std::vector<double> &sample);
    void FitUMVU(const std::vector<double> &sample);
};

#endif // UNIFORMRAND_H
