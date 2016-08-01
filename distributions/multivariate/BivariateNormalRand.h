#ifndef BIVARIATENORMALRAND_H
#define BIVARIATENORMALRAND_H

#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/NormalRand.h"

/**
 * @brief The BivariateNormalRand class
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
    std::string name() const override;

    void setLocation(DoublePair location);
    void setLocation(double location1, double location2);
    void setScale(const SquareMatrix<2> &rootCovariance);
    void setScale(double scale1, double scale2, double correlation);
    inline DoublePair getLocation() { return Mean(); }
    inline double getFirstLocation() const { return mu1; }
    inline double getSecondLocation() const { return mu2; }
    inline double getFirstScale() const { return sigma1; }
    inline double getSecondScale() const { return sigma2; }
    inline double getCorrelation() const { return rho; }

    double f(DoublePair point) const override;
    double F(DoublePair point) const override;
    DoublePair variate() const override;

    DoublePair Mean() const override;
    void Covariance(SquareMatrix<2> &matrix) const override;
    double Correlation() const override;

    void getFirstMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const;
    void getSecondMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const;
};

#endif // BIVARIATENORMALRAND_H
