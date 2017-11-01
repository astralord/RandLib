#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include "BetaRand.h"

/**
 * @brief The UniformRand class <BR>
 * Uniform continuous distribution
 *
 * f(x | a, b) = 1 / (b - a) for a < x < b
 *
 * Notation: X ~ U(a, b)
 *
 * Related distributions: <BR>
 * X ~ B(1, 1, a, b) <BR>
 * (X - a) / (b - a) ~ IH(1)
 */
class RANDLIBSHARED_EXPORT UniformRand : public BetaDistribution
{
public:
    UniformRand(double minValue = 0, double maxValue = 1);
    String Name() const override;

    using BetaDistribution::SetSupport;

    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

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

    static constexpr char TOO_LARGE_A[] = "Minimum element of the sample is smaller than lower boundary returned by method: ";
    static constexpr char TOO_SMALL_B[] = "Maximum element of the sample is larger than upper boundary returned by method: ";

public:
    inline double Entropy() const;

    /**
     * @fn FitMinimum
     * fit minimum with maximum-likelihood estimator if unbiased parameter is false,
     * otherwise fit minimum using UMVU estimator
     * @param sample
     */
    void FitMinimum(const std::vector<double> &sample, bool unbiased = false);
    /**
     * @fn FitMaximum
     * fit maximum with maximum-likelihood estimator if unbiased parameter is false,
     * otherwise fit maximum using UMVU estimator
     * @param sample
     */
    void FitMaximum(const std::vector<double> &sample, bool unbiased = false);
    /**
     * @fn Fit
     * fit support with maximum-likelihood estimator if unbiased parameter is false,
     * otherwise fit support using UMVU estimator
     * @param sample
     */
    void Fit(const std::vector<double> &sample, bool unbiased = false);
};

#endif // UNIFORMRAND_H
