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
    //TODO: find a way to initialize them without dummy
    /// Tables for ziggurat
    static double stairWidth[257], stairHeight[256];
    static double constexpr x1 = 7.69711747013104972;
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
    double F(double x) const override;
    double Variate() const override;

    static double Variate(double rate);
    double Median() const override;
    static double StandardVariate();

    std::complex<double> CF(double t) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double Entropy() const;

    double Moment(int n) const;

    /// Maximum-likelihood estimation
    bool FitMLE(const std::vector<double> &sample);
    /// Method of moments
    bool FitMM(const std::vector<double> &sample);
    /// Uniformly minimum variance unbiased estimator
    bool FitUMVU(const std::vector<double> &sample);
};

#endif // EXPONENTIALRAND_H
