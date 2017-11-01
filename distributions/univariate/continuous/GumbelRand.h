#ifndef GUMBELRAND_H
#define GUMBELRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GumbelRand class <BR>
 * Gumbel distribution
 *
 * Notation: X ~ Gumbel(μ, β)
 *
 * Related distributions: <BR>
 * exp(-(X - μ) / β) ~ Exp(1)
 */
class RANDLIBSHARED_EXPORT GumbelRand : public ContinuousDistribution
{
    double mu = 0; ///< location μ
    double beta = 1; ///< scale β
    double logBeta = 0; ///< log(β)
public:
    GumbelRand(double location, double scale);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return beta; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
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

public:
    double Entropy() const;
};

#endif // GUMBELRAND_H
