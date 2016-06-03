#ifndef NORMALINVERSEGAMMARAND_H
#define NORMALINVERSEGAMMARAND_H

#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/InverseGammaRand.h"

/**
 * @brief The NormalInverseGammaRand class
 */
class RANDLIBSHARED_EXPORT NormalInverseGammaRand : public BivariateProbabilityDistribution<double, double>
{
    double mu, lambda;
    InverseGammaRand Y;
    double pdfCoef, cdfCoef;

public:
    NormalInverseGammaRand(double location = 0, double precision = 1, double shape = 1, double rate = 1);
    std::string name() override;

    void setParameters(double location, double precision, double shape, double rate);
    inline double getLocation() { return mu; }
    inline double getPrecision() { return lambda; }
    inline double getShape() { return Y.getShape(); }
    inline double getRate() { return Y.getRate(); }

    double f(DoublePair point) const override;
    double F(DoublePair point) const override;
    DoublePair variate() const override;

    DoublePair Mean() const override;
    void Covariance(SquareMatrix<2> &matrix) const override;
    double Correlation() const override;

    void getFirstMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const;
    void getSecondMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const;
};

#endif // NORMALINVERSEGAMMARAND_H
