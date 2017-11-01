#ifndef BIVARIATENORMALRAND_H
#define BIVARIATENORMALRAND_H

#include "ContinuousBivariateDistribution.h"
#include "../univariate/continuous/NormalRand.h"

/**
 * @brief The BivariateNormalRand class <BR>
 * Bivariate Gaussian (normal) distribution
 *
 * Notation: X ~ N(μ1, μ2, σ1, ρ, σ2)
 */
class RANDLIBSHARED_EXPORT BivariateNormalRand : public ContinuousBivariateDistribution<NormalRand, NormalRand>
{
    double mu1 = 0; ///< first location μ1
    double mu2 = 0; ///< second location μ2
    double sigma1 = 1; ///< first scale σ1
    double sigma2 = 1; ///< second scale σ2
    double rho = 0; ///< correlation coefficient ρ
    double pdfCoef = 1 + M_LN2 + M_LNPI; ///< coefficient for faster pdf calculation
    double sqrt1mroSq = 1; ///< √(1 - ρ)

public:
    BivariateNormalRand(double location1, double location2, double scale1, double scale2, double correlation);
    String Name() const override;

    void SetLocations(double location1, double location2);
    void SetCovariance(double scale1, double scale2, double correlation);
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
