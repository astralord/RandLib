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
    NormalRand X;
    UniformRand U;

    double pdfCoef; /// sqrt(l / (2 * pi))
    double cdfCoef; /// exp(2 * l / mu)
public:
    WaldRand(double mean = 1, double shape = 1);

    void setParameters(double mean, double shape);
    inline double getMean() const { return mu; }
    inline double getShape() const { return l; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return mu; }
    double Var() const override { return mu * mu * mu / l; }
};


#endif // WALDRAND_H
