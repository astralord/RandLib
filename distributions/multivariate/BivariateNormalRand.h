#ifndef BIVARIATENORMALRAND_H
#define BIVARIATENORMALRAND_H

#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/NormalRand.h"

/**
 * @brief The BivariateNormalRand class
 */
class RANDLIBSHARED_EXPORT BivariateNormalRand : public BivariateProbabilityDistribution
{
    double mu1, mu2;
    double sigma1, sigma2, ro;

    NormalRand X, Y;
    double pdfCoef, sqrt1mroSq;

public:
    BivariateNormalRand(DoublePair location, SquareMatrix<2> rootCovariance);
    BivariateNormalRand(double mu1, double mu2, double sigma1, double sigma2, double correlation);
    std::string name() override;

    void setLocation(DoublePair location);
    void setLocation(double location1, double location2);
    void setScale(SquareMatrix<2> rootCovariance);
    void setScale(double scale1, double scale2, double correlation);
    inline DoublePair getLocation() { return Mean(); }
    inline double getFirstLocation() { return mu1; }
    inline double getSecondLocation() { return mu2; }
    inline double getFirstScale() { return sigma1; }
    inline double getSecondScale() { return sigma2; }
    inline double getCorrelation() { return ro; }

    double f(DoublePair point) const override;
    double F(DoublePair point) const override;
    DoublePair variate() const override;

    DoublePair Mean() const override;
    void Covariance(SquareMatrix<2> &matrix) const override;
    double Correlation() const override;

    void getFirstConditionalDistribution(NormalRand & distribution, double y);
    void getSecondConditionalDistribution(NormalRand & distribution, double x);

    void getFirstMarginalDistribution(ProbabilityDistribution<double> &distribution) const;
    void getSecondMarginalDistribution(ProbabilityDistribution<double> &distribution) const;
};

#endif // BIVARIATENORMALRAND_H
