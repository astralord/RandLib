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
class RANDLIBSHARED_EXPORT BivariateNormalRand : public BivariateProbabilityDistribution
{
    double mu1, mu2;
    double sigma1, sigma2, rho;

    NormalRand X, Y;
    double pdfCoef, sqrt1mroSq;

public:
    BivariateNormalRand(const DoublePair &location, const DoubleTriplet &covariance);
    BivariateNormalRand(double location1, double location2, double scale1, double correlation, double scale2);
    std::string Name() const override;
    DoublePair MinValue() const { return DoublePair(-INFINITY, -INFINITY); }
    DoublePair MaxValue() const { return DoublePair(INFINITY, INFINITY); }

    void SetLocation(const DoublePair &location);
    void SetLocation(double location1, double location2);
    void SetScale(const DoubleTriplet &covariance);
    void SetScale(double scale1, double correlation, double scale2);
    inline DoublePair GetLocation() { return Mean(); }
    inline double GetFirstLocation() const { return mu1; }
    inline double GetSecondLocation() const { return mu2; }
    inline double GetFirstScale() const { return sigma1; }
    inline double GetSecondScale() const { return sigma2; }
    inline double GetCorrelation() const { return rho; }

    double f(const DoublePair &point) const override;
    double F(const DoublePair & point) const override;
    DoublePair Variate() const override;

    DoublePair Mean() const override;
    DoubleTriplet Covariance() const override;
    double Correlation() const override;

    void GetMarginalDistributions(ContinuousDistribution &distribution1, ContinuousDistribution &distribution2) const override;
};

#endif // BIVARIATENORMALRAND_H
