#ifndef WRAPPEDEXPONENTIALRAND_H
#define WRAPPEDEXPONENTIALRAND_H

#include "CircularDistribution.h"

/**
 * @brief The WrappedExponentialRand class <BR>
 * Wrapped Exponential distribution
 *
 * Notation: X ~ WE(λ)
 */
class RANDLIBSHARED_EXPORT WrappedExponentialRand : public CircularDistribution
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

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;

    double CircularMean() const override;
    double CircularVariance() const override;
    double Median() const override;
    double Mode() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // WRAPPEDEXPONENTIALRAND_H
