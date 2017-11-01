#ifndef DEGENERATERAND_H
#define DEGENERATERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The DegenerateRand class <BR>
 * Degenerate distribution
 *
 * f(x|a) = δ(a)
 *
 * Notation: X ~ δ(a)
 */
class RANDLIBSHARED_EXPORT DegenerateRand : public ContinuousDistribution
{
    double a = 0; ///< value

public:
    explicit DegenerateRand(double value = 0);
    String Name() const override;
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
     * @fn Fit
     * @param sample
     */
    void Fit(const std::vector<double> &sample);
};

#endif // DEGENERATERAND_H
