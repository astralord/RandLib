#ifndef DEGENERATERAND_H
#define DEGENERATERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The DegenerateRand class
 * Degenerate distribution
 *
 * f(x|a) = δ(a)
 *
 * Notation: X ~ δ(a)
 */
class RANDLIBSHARED_EXPORT DegenerateRand : public ContinuousDistribution
{
    double a;

public:
    explicit DegenerateRand(double value = 0);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return a; }

    void SetValue(double value);
    inline double GetValue() const { return a; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double Variate() const override;

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

    double Entropy() const;

    /**
     * @brief FitMLE
     * @param sample
     * @return
     */
    void FitMLE(const std::vector<double> &sample);
    /**
     * @brief FitMM
     * @param sample
     * @return
     */
    void FitMM(const std::vector<double> &sample);
};

#endif // DEGENERATERAND_H
