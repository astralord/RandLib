#ifndef NORMALINVERSEGAMMARAND_H
#define NORMALINVERSEGAMMARAND_H

#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/NormalRand.h"
#include "../univariate/continuous/InverseGammaRand.h"

/**
 * @brief The NormalInverseGammaRand class
 */
class RANDLIBSHARED_EXPORT NormalInverseGammaRand : public BivariateProbabilityDistribution
{
    double mu, lambda;
    InverseGammaRand Y;
    double pdfCoef;

public:
    NormalInverseGammaRand(double location = 0, double precision = 1, double shape = 1, double rate = 1);
    std::string name();

    void setParameters(double location, double precision, double shape, double rate);
    inline double getLocation() { return mu; }
    inline double getPrecision() { return lambda; }
    inline double getShape() { return Y.getShape(); }
    inline double getRate() { return Y.getRate(); }

    double f(double2d grid) const;
    double F(double2d grid) const;
    double2d variate() const;

    double2d Mean() const;
};

#endif // NORMALINVERSEGAMMARAND_H
