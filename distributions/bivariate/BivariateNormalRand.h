#ifndef BIVARIATENORMALRAND_H
#define BIVARIATENORMALRAND_H

#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/NormalRand.h"

/**
 * @brief The BivariateNormalRand class
 * Bivariate Gaussian (normal) distribution
 *
 * Notation: X ~ N(μ_1, μ_2, σ_1, ρ, σ_2)
 */
class RANDLIBSHARED_EXPORT BivariateNormalRand : public BivariateProbabilityDistribution<NormalRand, NormalRand>
{
    double mu1, mu2;
    double sigma1, sigma2, rho;
    double pdfCoef, sqrt1mroSq;

public:
    BivariateNormalRand(double location1, double location2, double scale1, double correlation, double scale2);
    std::string Name() const override;
    DoublePair MinValue() const { return DoublePair(-INFINITY, -INFINITY); }
    DoublePair MaxValue() const { return DoublePair(INFINITY, INFINITY); }

    void SetLocations(double location1, double location2);
    void SetCovariance(double scale1, double correlation, double scale2);
    inline DoublePair GetLocation() { return Mean(); }
    inline double GetFirstLocation() const { return mu1; }
    inline double GetSecondLocation() const { return mu2; }
    inline double GetFirstScale() const { return sigma1; }
    inline double GetSecondScale() const { return sigma2; }
    inline double GetCorrelation() const { return rho; }

    double f(const DoublePair &point) const override;
    double logf(const DoublePair &point) const override;
    double F(const DoublePair & point) const override;
    DoublePair Variate() const override;

    double Correlation() const override;
};

#endif // BIVARIATENORMALRAND_H
