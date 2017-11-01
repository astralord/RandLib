#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "GammaRand.h"

/**
 * @brief The ExponentialRand class <BR>
 * Exponential distribution
 *
 * f(x | β) = β exp(-βx)
 *
 * Notation: X ~ Exp(β)
 *
 * Related distributions: <BR>
 * X ~ Γ(1, β)
 */
class RANDLIBSHARED_EXPORT ExponentialRand : public FreeScaleGammaDistribution
{
    /// Tables for ziggurat
    static long double stairWidth[257], stairHeight[256];
    static constexpr long double x1 = 7.69711747013104972l;
    static bool dummy;
    static bool SetupTables();

public:
    explicit ExponentialRand(double rate = 1) : FreeScaleGammaDistribution(1, rate) {}

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;
    static double StandardVariate();

    double Median() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;
    double Moment(int n) const;
    double ThirdMoment() const override { return Moment(3); }
    double FourthMoment() const override { return Moment(4); }
};

#endif // EXPONENTIALRAND_H
