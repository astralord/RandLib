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
 * X ~ Beta(1, 1, a, b) <BR>
 * (X - a) / (b - a) ~ IH(1)
 */
class RANDLIBSHARED_EXPORT UniformRand : public BetaDistribution
{
public:
    UniformRand(double minValue = 0, double maxValue = 1);
    std::string Name() const override;

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
     * @fn FitMinimumMLE
     * fit minimum with maximum-likelihood estimator
     * @param sample
     */
    void FitMinimumMLE(const std::vector<double> &sample);
    /**
     * @fn FitMaximumMLE
     * fit maximum with maximum-likelihood estimator
     * @param sample
     */
    void FitMaximumMLE(const std::vector<double> &sample);
    /**
     * @fn FitSupportMLE
     * fit support with maximum-likelihood estimator
     * @param sample
     */
    void FitSupportMLE(const std::vector<double> &sample);

    /**
     * @fn FitMinimumMM
     * fit minimum with method of moments
     * @param sample
     */
    void FitMinimumMM(const std::vector<double> &sample);
    /**
     * @fn FitMaximumMM
     * fit maximum with method of moments
     * @param sample
     */
    void FitMaximumMM(const std::vector<double> &sample);
    /**
     * @fn FitSupportMM
     * fit support with method of moments
     * @param sample
     */
    void FitSupportMM(const std::vector<double> &sample);
    
    /**
     * @fn FitMinimumUMVU
     * fit minimum via uniformly-minimum variance unbiased estimator
     * @param sample
     */
    void FitMinimumUMVU(const std::vector<double> &sample);
    /**
     * @fn FitMaximumUMVU
     * fit maximum via uniformly-minimum variance unbiased estimator
     * @param sample
     */
    void FitMaximumUMVU(const std::vector<double> &sample);
    /**
     * @fn FitSupportUMVU
     * fit support via uniformly-minimum variance unbiased estimator
     * @param sample
     */
    void FitSupportUMVU(const std::vector<double> &sample);
};

#endif // UNIFORMRAND_H
