#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "GammaRand.h"

/**
 * @brief The ExponentialRand class
 * Exponential distribution
 * X ~ Exp(\lambda)
 *
 * f(x|\lambda) = \lambda \exp(-\lambda x)
 *
 * X ~ \Gamma(1, 1 / \lambda)
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
    inline double getRate() const { return beta; }

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

    bool fit_MM(const QVector<double> &sample) override;
    bool fit_MLE(const QVector<double> &sample) override;
    bool fit_UMVU(const QVector<double> &sample);
};

#endif // EXPONENTIALRAND_H
