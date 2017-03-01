#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "GammaRand.h"

/**
 * @brief The ExponentialRand class
 * Exponential distribution
 *
 * f(x | β) = β exp(-βx)
 *
 * Notation: X ~ Exp(β)
 *
 * Related distributions:
 * X ~ Gamma(1, β)
 */
class RANDLIBSHARED_EXPORT ExponentialRand : public GammaRand
{
    /// Tables for ziggurat
    static long double stairWidth[257], stairHeight[256];
    static constexpr long double x1 = 7.69711747013104972l;
    static bool dummy;
    static bool SetupTables();

public:
    explicit ExponentialRand(double rate = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

protected:
    using GammaRand::SetParameters;

public:
    void SetRate(double rate);

    double f(double x) const override;
    double logf(double x) const override;
    double F(double x) const override;
    double S(double x) const override;
    double Variate() const override;
    static double StandardVariate();

    double Median() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;

    double Moment(int n) const;

    /// Maximum-likelihood estimation
    void FitMLE(const std::vector<double> &sample);
    /// Method of moments
    void FitMM(const std::vector<double> &sample);
    /// Uniformly minimum variance unbiased estimator
    void FitUMVU(const std::vector<double> &sample);
};

#endif // EXPONENTIALRAND_H
