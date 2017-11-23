#ifndef NONCENTRALCHISQUAREDRAND_H
#define NONCENTRALCHISQUAREDRAND_H

#include "ContinuousDistribution.h"
#include "GammaRand.h"
#include "NormalRand.h"
#include "../discrete/PoissonRand.h"

/**
 * @brief The NoncentralChiSquaredRand class <BR>
 * Noncentral Chi-Squared distribution
 *
 * Notation: X ~ χ'^2(k, λ)
 *
 * Related distributions: <BR>
 * If X ~ χ'^2(k, 0), then X ~ χ^2(k) <BR>
 * X ~ χ^2(k + 2J), where J ~ Po(λ)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT NoncentralChiSquaredRand : public ContinuousDistribution<RealType>
{
    double k = 1; ///< degree
    double lambda = 2; ///< noncentrality λ
    double halfK = 0.5; ///< k / 2
    double halfLambda = 1; ///< λ / 2
    double sqrtLambda = M_SQRT2; ///< √λ
    double logLambda = M_LN2; ///< log(λ)

    PoissonRand<int> Y{};

public:
    explicit NoncentralChiSquaredRand(double degree = 1, double noncentrality = 0);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double degree, double noncentrality);
    inline double GetDegree() const { return k; }
    inline double GetNoncentrality() const { return lambda; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

private:
    RealType variateForDegreeEqualOne() const;

public:
    static RealType Variate(double degree, double noncentrality, RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // NONCENTRALCHISQUAREDRAND_H
