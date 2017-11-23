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
template < typename RealType = double >
class RANDLIBSHARED_EXPORT LogisticRand : public ContinuousDistribution<RealType>
{
    double mu = 0; ///< location μ
    double s = 1; ///< scale s
    double logS = 0; ///< log(s)

public:
    LogisticRand(double location = 0, double scale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    RealType MinValue() const override { return -INFINITY; }
    RealType MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return s; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;

    /**
     * @fn FitLocation
     * fit location parameter via maximum-likelihood
     * @param sample
     */
    void FitLocation(const std::vector<RealType> &sample);
};

#endif // LOGISTICRAND_H
