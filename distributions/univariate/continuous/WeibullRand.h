#ifndef WEIBULLRAND_H
#define WEIBULLRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The WeibullRand class <BR>
 * Weibull distribution
 *
 * Notation: X ~ Weibull(λ, k)
 */
class RANDLIBSHARED_EXPORT WeibullRand : public ContinuousDistribution
{
    double lambda = 1; ///< scale λ
    double k = 1; ///< shape k
    double kInv = 1; /// 1 /k
    double logk_lambda = 0; /// log(k/λ)

public:
    WeibullRand(double scale = 1, double shape = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double scale, double shape);
    inline double GetScale() const { return lambda; }
    inline double GetShape() const { return k; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
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
};

#endif // WEIBULLRAND_H
