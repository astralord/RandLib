#ifndef WALDRAND_H
#define WALDRAND_H

#include "ContinuousRand.h"
#include "NormalRand.h"
#include "UniformRand.h"

/**
 * @brief The WaldRand class
 * Wald (Inverse Gaussian) distribution
 * X ~ IG(mu, l)
 */
class RANDLIBSHARED_EXPORT WaldRand : public ContinuousRand
{
    double mu, l;

    double pdfCoef; /// sqrt(l / (2 * pi))
    double cdfCoef; /// exp(2 * l / mu)
public:
    WaldRand(double mean = 1, double shape = 1);
    std::string name() override;

    void setParameters(double mean, double shape);
    inline double getMean() const { return mu; }
    inline double getShape() const { return l; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override { return mu; }
    double Variance() const override { return mu * mu * mu / l; }
    
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


#endif // WALDRAND_H
