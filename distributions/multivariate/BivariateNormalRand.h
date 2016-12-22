#ifndef BIVARIATENORMALRAND_H
#define BIVARIATENORMALRAND_H

#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/NormalRand.h"

/**
 * @brief The BivariateNormalRand class
 * Bivariate Gaussian (normal) distribution
 *
 * Notation: X ~ N(μ_1, μ_2, σ_1, σ_2, ρ)
 */
class RANDLIBSHARED_EXPORT BivariateNormalRand : public BivariateProbabilityDistribution<double, double>
{
    double mu1, mu2;
    double sigma1, sigma2, rho;

    NormalRand X, Y;
    double pdfCoef, sqrt1mroSq;

public:
    BivariateNormalRand(DoublePair location, const SquareMatrix<2> &rootCovariance);
    BivariateNormalRand(double mu1, double mu2, double sigma1, double sigma2, double correlation);
    std::string Name() const override;
    DoublePair MinValue() const { return DoublePair(-INFINITY, -INFINITY); }
    DoublePair MaxValue() const { return DoublePair(INFINITY, INFINITY); }

    void SetLocation(DoublePair location);
    void SetLocation(double location1, double location2);
    void SetScale(const SquareMatrix<2> &rootCovariance);
    void SetScale(double scale1, double scale2, double correlation);
    inline DoublePair GetLocation() { return Mean(); }
    inline double GetFirstLocation() const { return mu1; }
    inline double GetSecondLocation() const { return mu2; }
    inline double GetFirstScale() const { return sigma1; }
    inline double GetSecondScale() const { return sigma2; }
    inline double GetCorrelation() const { return rho; }

    double f(DoublePair point) const override;
    double F(DoublePair point) const override;
    DoublePair Variate() const override;

    DoublePair Mean() const override;
    void Covariance(SquareMatrix<2> &matrix) const override;
    double Correlation() const override;

    void GetFirstMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const;
    void GetSecondMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const;
};

#endif // BIVARIATENORMALRAND_H
