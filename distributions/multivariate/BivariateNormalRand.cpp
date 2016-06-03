#include "BivariateNormalRand.h"


BivariateNormalRand::BivariateNormalRand(DoublePair location, SquareMatrix<2> rootCovariance)
{
    setLocation(location);
    setScale(rootCovariance);
}

BivariateNormalRand::BivariateNormalRand(double mu1, double mu2, double sigma1, double sigma2, double correlation)
{
    setLocation(mu1, mu2);
    setScale(sigma1, sigma2, correlation);
}

std::string BivariateNormalRand::name()
{
    return "Bivariate Normal( (" + toStringWithPrecision(getFirstLocation()) + ", "
                               + toStringWithPrecision(getSecondLocation()) + "), ("
                               + toStringWithPrecision(getFirstScale()) + ", "
                               + toStringWithPrecision(getSecondScale()) + ", "
                               + toStringWithPrecision(getCorrelation()) + ") )";
}

void BivariateNormalRand::setLocation(DoublePair location)
{
    mu1 = location.first;
    mu2 = location.second;
    X.setLocation(mu1);
    Y.setLocation(mu2);
}

void BivariateNormalRand::setLocation(double location1, double location2)
{
    mu1 = location1;
    mu2 = location2;
    X.setLocation(mu1);
    Y.setLocation(mu2);
}

void BivariateNormalRand::setScale(double scale1, double scale2, double correlation)
{
    sigma1 = scale1;
    if (sigma1 <= 0)
        sigma1 = 1.0;

    sigma2 = scale2;
    if (sigma2 <= 0)
        sigma2 = 1.0;

    X.setScale(sigma1);
    Y.setScale(sigma2);

    ro = correlation;
    if (ro < 0 || ro >= 1) // we don't accept ro = 1 for now (because of pdf)
        ro = 0.0;

    sqrt1mroSq = std::sqrt(1.0 - ro * ro);
    pdfCoef = 0.5 * M_1_PI / (sigma1 * sigma2 * sqrt1mroSq);
}

void BivariateNormalRand::setScale(SquareMatrix<2> rootCovariance)
{
    setScale(rootCovariance(0, 0), rootCovariance(1, 1), rootCovariance(1, 0));
}

double BivariateNormalRand::f(DoublePair point) const
{
    if (ro == 0.0)
        return X.f(point.first) * Y.f(point.second);
    if (ro == 1.0)
        return NAN; // TODO!!!
    double xAdj = (point.first - mu1) / sigma1, yAdj = (point.second - mu2) / sigma2;
    double z = xAdj * xAdj - 2 * ro * xAdj * yAdj + yAdj * yAdj;
    return pdfCoef * std::exp(-0.5 * z / (1 - ro * ro));
}

double BivariateNormalRand::F(DoublePair point) const
{
    if (ro == 0.0)
    {
        return X.F(point.first) * Y.F(point.second);
    }
    if (std::fabs(ro) == 1.0)
    {
        double secondPoint = ro * (point.second - mu2) * sigma1 / sigma2 + mu1;
        return X.F(std::min(point.first, secondPoint));
    }
    // TODO!!
    return NAN;
}

DoublePair BivariateNormalRand::variate() const
{
    double Z1 = NormalRand::standardVariate();
    double Z2 = NormalRand::standardVariate();
    double x = mu1 + sigma1 * Z1;
    double y = mu2 + sigma2 * (ro * Z1 + sqrt1mroSq * Z2);
    return std::make_pair(x, y);
}

DoublePair BivariateNormalRand::Mean() const
{
    return std::make_pair(mu1, mu2);
}

void BivariateNormalRand::Covariance(SquareMatrix<2> &matrix) const
{
    matrix(0, 0) = sigma1 * sigma1;
    matrix(0, 1) = matrix(1, 0) = ro * sigma1 * sigma2;
    matrix(1, 1) = sigma2 * sigma2;
}

double BivariateNormalRand::Correlation() const
{
    return ro;
}

void BivariateNormalRand::getFirstMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const
{
    distribution = X;
}

void BivariateNormalRand::getSecondMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const
{
    distribution = Y;
}

