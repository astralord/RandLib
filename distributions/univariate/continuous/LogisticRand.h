#ifndef LOGISTICRAND_H
#define LOGISTICRAND_H

#include "ContinuousDistribution.h"
#include "UniformRand.h"

/**
 * @brief The LogisticRand class <BR>
 * Logistic distribution
 *
 * Notation: X ~ Logistic(μ, s)
 *
 * Related distributions: <BR>
 * 1 / (exp((X - μ) / s) + 1) ~ U(0, 1)
 */
class RANDLIBSHARED_EXPORT LogisticRand : public ContinuousDistribution
{
    double mu = 0; ///< location μ
    double s = 1; ///< scale s
    double logS = 0; ///< log(s)

public:
    LogisticRand(double location = 0, double scale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return s; }

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

    /**
     * @brief FitLocation
     * fit location parameter via maximum-likelihood
     * @param sample
     */
    void FitLocation(const std::vector<double> &sample);
};

#endif // LOGISTICRAND_H
