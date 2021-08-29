#ifndef NORMALINVERSEGAMMARAND_H
#define NORMALINVERSEGAMMARAND_H

#include "ContinuousBivariateDistribution.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/InverseGammaRand.h"

/**
 * @brief The NormalInverseGammaRand class <BR>
 * Normal-inverse-gamma distribution
 *
 * Notation: X ~ NIG(μ, λ, α, β)
 */
template< typename RealType = double>
class RANDLIBSHARED_EXPORT NormalInverseGammaRand : public ContinuousBivariateDistribution<StudentTRand<RealType>, InverseGammaRand<RealType>, RealType>
{
    double mu = 0; ///< location μ
    double lambda = 1; ///< precision λ
    double alpha = 1; ///< shape α
    double beta = 1; ///< rate β
    double pdfCoef = 0.5 * (M_LNPI - M_LN2); ///< coefficient for faster pdf calculation

public:
    NormalInverseGammaRand(double location = 0, double precision = 1, double shape = 1, double rate = 1);
    String Name() const override;

    void SetParameters(double location, double precision, double shape, double rate);
    inline double GetLocation() const { return mu; }
    inline double GetPrecision() const { return lambda; }
    inline double GetShape() const { return alpha; }
    inline double GetRate() const { return beta; }

    double f(const Pair<RealType> &point) const override;
    double logf(const Pair<RealType> &point) const override;
    double F(const Pair<RealType> & point) const override;
    Pair<RealType> Variate() const override;

    long double Correlation() const override;

    Pair<RealType> Mode() const override;
};

#endif // NORMALINVERSEGAMMARAND_H
