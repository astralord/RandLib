#ifndef WEIBULLRAND_H
#define WEIBULLRAND_H

#include "ContinuousDistribution.h"
#include "InverseGammaRand.h"

/**
 * @brief The WeibullRand class <BR>
 * Weibull distribution
 *
 * Notation: X ~ Weibull(λ, k)
 */
template < typename RealType = long double >
class RANDLIBSHARED_EXPORT WeibullRand : public ContinuousDistribution<RealType>
{
    double lambda = 1; ///< scale λ
    double k = 1; ///< shape k
    double kInv = 1; /// 1 /k
    double logk_lambda = 0; /// log(k/λ)

public:
    WeibullRand(double scale = 1, double shape = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double scale, double shape);
    inline double GetScale() const { return lambda; }
    inline double GetShape() const { return k; }

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

protected:
    double getPowSampleMean(const std::vector<RealType> &sample) const;

    /**
     * @fn FitScale
     * Fit λ by maximum-likelihood
     * @param sample
     */
    void FitScale(const std::vector<RealType> &sample);

    /**
     * @fn FitScaleBayes
     * Fit λ, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    InverseGammaRand<RealType> FitScaleBayes(const std::vector<RealType> &sample, const InverseGammaRand<RealType> &priorDistribution, bool MAP = false);
};

#endif // WEIBULLRAND_H
