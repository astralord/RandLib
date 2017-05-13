#ifndef KUMARASWAMYRAND_H
#define KUMARASWAMYRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The KumaraswamyRand class <BR>
 * Kumaraswamy distribution
 *
 * f(x|a,b) = abx^{a-1} (1-x^a)^{b-1}
 *
 * Notation: X ~ Kumaraswamy(a, b)
 */
class RANDLIBSHARED_EXPORT KumaraswamyRand : public ContinuousDistribution
{
    double a = 1; ///< first shape
    double b = 1; ///< second shape
    double logA = 0; ///< log(a)
    double logB = 0; ///< log(b)

public:
    explicit KumaraswamyRand(double shape1, double shape2);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return 1; }

    void SetShapes(double shape1, double shape2);
    inline double GetFirstShape() const { return a; }
    inline double GetSecondShape() const { return b; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Moment(int n) const;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};

#endif // KUMARASWAMYRAND_H
