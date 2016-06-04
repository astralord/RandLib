#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "GammaRand.h"

/**
 * @brief The ExponentialRand class
 * Exponential distribution
 * X ~ Exp(\beta)
 *
 * f(x|\beta) = \beta \exp(-\beta x)
 *
 * X ~ \Gamma(1, 1 / \beta)
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
    std::string name() override;

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

    double Quantile(double p) const override;

    double Median() const override;

    double Entropy() const;

    double Moment(int n) const;

    bool fitRateMM(const std::vector<double> &sample);
    bool fitRateMLE(const std::vector<double> &sample);
    bool fitRateUMVU(const std::vector<double> &sample);
};

#endif // EXPONENTIALRAND_H
