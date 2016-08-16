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
    static bool setupTables();

public:
    explicit ExponentialRand(double rate = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

protected:
    using GammaRand::setParameters;

public:
    void setRate(double rate);

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double rate);
    static double standardVariate();

    std::complex<double> CF(double t) const override;

    double QuantileImpl(double p) const override;

    double Median() const override;

    double Entropy() const;

    double Moment(int n) const;

    /// Maximum-likelihood estimation
    bool fitMLE(const std::vector<double> &sample);
    /// Method of moments
    bool fitMM(const std::vector<double> &sample);
    /// Uniformly minimum variance unbiased estimator
    bool fitUMVU(const std::vector<double> &sample);
};

#endif // EXPONENTIALRAND_H
