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

    /// Delete GammaRand estimators for shape
    void FitShapeAndRateMLE(const std::vector<double> &sample) = delete;
    void FitShapeMM(const std::vector<double> &sample) = delete;
    void FitShapeAndRateMM(const std::vector<double> &sample) = delete;
};

#endif // EXPONENTIALRAND_H
