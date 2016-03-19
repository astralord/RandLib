#ifndef LEVYRAND_H
#define LEVYRAND_H

#include "StableRand.h"

/**
 * @brief The LevyRand class
 *
 * f(x|\mu, c) = \sqrt(c exp(c / (\mu - x)) / (2 \pi (x - \mu)^3)
 *
 * X ~ Levy(\mu, c)
 * X ~ Stable(0.5, 1, \mu, c)
 */
class RANDLIBSHARED_EXPORT LevyRand : public StableRand
{
    double c;
    double pdfCoef; /// sqrt(0.5 * sigma / pi)

public:
    LevyRand(double location = 0, double scale = 1);
    std::string name() override;

protected:
    void setParameters(double exponent, double skewness, double scale, double location);

public:
    void setParameters(double location, double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return sigma; }

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
