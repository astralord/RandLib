#ifndef WRAPPEDEXPONENTIALRAND_H
#define WRAPPEDEXPONENTIALRAND_H

#include "CircularDistribution.h"

/**
 * @brief The WrappedExponentialRand class <BR>
 * Wrapped Exponential distribution
 *
 * Notation: X ~ WE(λ)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT WrappedExponentialRand : public CircularDistribution<RealType>
{
    double lambda = 1; ///< rate λ
    double logLambda = 0; ///< ln(λ)
    double scaledLambda = 2 * M_PI; ///< 2πλ
    double pdfCoef = 0, logpdfCoef = 0, expmScaledLambda = 0;

public:
    WrappedExponentialRand(double rate);

    String Name() const override;

    void SetRate(double rate);
    inline double GetRate() const { return lambda; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    long double CircularMean() const override;
    long double CircularVariance() const override;
    RealType Median() const override;
    RealType Mode() const override;

protected:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // WRAPPEDEXPONENTIALRAND_H
