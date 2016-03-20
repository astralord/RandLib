#ifndef LEVYRAND_H
#define LEVYRAND_H

#include "StableRand.h"

/**
 * @brief The LevyRand class
 * Levy distribution
 * X ~ Levy(\mu, \sigma)
 *
 * f(x|\mu, \sigma) = \sqrt(\sigma exp(\sigma / (\mu - x)) / (2 \pi (x - \mu)^3)
 *
 * X ~ Stable(0.5, 1, \sigma, \mu)
 */
class RANDLIBSHARED_EXPORT LevyRand : public StableRand
{
public:
    LevyRand(double location = 0, double scale = 1);
    std::string name() override;

private:
    using StableRand::setParameters;

public:
    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double scale);
    static double standardVariate();
    
    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;
    
    double Mode() const override;
    
    /// Verify that all elements of sample can have this distribution
    bool checkValidity(const QVector<double> &sample);
    
    /// Maximum likelihood estimators
    bool fitScale_MLE(const QVector<double> &sample);
};

#endif // LEVYRAND_H
