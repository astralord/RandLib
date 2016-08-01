#include "BivariateNormalRand.h"


BivariateNormalRand::BivariateNormalRand(DoublePair location, const SquareMatrix<2> &rootCovariance)
{
    setLocation(location);
    setScale(rootCovariance);
}

BivariateNormalRand::BivariateNormalRand(double mu1, double mu2, double sigma1, double sigma2, double correlation)
{
    setLocation(mu1, mu2);
    setScale(sigma1, sigma2, correlation);
}

std::string BivariateNormalRand::name() const
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

    rho = correlation;
    if (rho < 0 || rho >= 1) // we don't accept ro = 1 for now (because of pdf)
        rho = 0.0;

    sqrt1mroSq = std::sqrt(1.0 - rho * rho);
    pdfCoef = 0.5 * M_1_PI / (sigma1 * sigma2 * sqrt1mroSq);
}

void BivariateNormalRand::setScale(const SquareMatrix<2> &rootCovariance)
{
    setScale(rootCovariance(0, 0), rootCovariance(1, 1), rootCovariance(1, 0));
}

double BivariateNormalRand::f(DoublePair point) const
{
    if (rho == 0.0) /// f(x, y) = f(x)f(y)
        return X.f(point.first) * Y.f(point.second);
    double xAdj = (point.first - mu1) / sigma1, yAdj = (point.second - mu2) / sigma2;
    if (rho == 1.0) /// f(x, y) = δ(xAdj - yAdj)
        return (xAdj - yAdj == 0) ? INFINITY : 0.0;
    if (rho == -1.0) /// f(x, y) = δ(xAdj + yAdj)
        return (xAdj + yAdj == 0) ? INFINITY : 0.0;
    double z = xAdj * xAdj - 2 * rho * xAdj * yAdj + yAdj * yAdj;
    return pdfCoef * std::exp(-0.5 * z / (1 - rho * rho));
}

double BivariateNormalRand::F(DoublePair point) const
{
    if (rho == 0.0) /// F(x, y) = F(x)F(y)
        return X.F(point.first) * Y.F(point.second);
    if (std::fabs(rho) == 1.0) {
        double secondPoint = rho * (point.second - mu2) * sigma1 / sigma2 + mu1;
        return X.F(std::min(point.first, secondPoint));
    }

    /// Unnumbered equation between (3) and (4) in Section 2.2 of Genz (2004),
    /// integrating in terms of theta between asin(ρ) and +/- π/2
    double p1 = 0;
    double x = point.first, y = point.second;
    double xAdj = (x - mu1) / sigma1, yAdj = (y - mu2) / sigma2;
    if (rho > 0) {
        p1 = (xAdj < yAdj) ? X.F(xAdj) : Y.F(yAdj);
    }
    else {
        p1 = std::max(X.F(xAdj) - Y.F(-yAdj), 0.0);
    }
    double lowLimit = std::asin(rho);
    double highLimit = RandMath::sign(rho) * M_PI_2;
    double p2 = RandMath::integral([this, x, y] (double theta) {
        /// Integrand is exp(-(x^2 + y^2 - 2xysin(θ)) / (2cos(θ)^2))
        double sinTheta = std::sin(theta);
        double cosTheta = std::cos(theta);
        double integrand = (x * sinTheta - y);
        integrand *= integrand;
        integrand /= (cosTheta * cosTheta);
        return std::exp(-0.5 * integrand);
    }, lowLimit, highLimit);
    return p1 - 0.5 * p2 / M_PI;
}

DoublePair BivariateNormalRand::variate() const
{
    double Z1 = NormalRand::standardVariate();
    double Z2 = NormalRand::standardVariate();
    double x = mu1 + sigma1 * Z1;
    double y = mu2 + sigma2 * (rho * Z1 + sqrt1mroSq * Z2);
    return std::make_pair(x, y);
}

DoublePair BivariateNormalRand::Mean() const
{
    return std::make_pair(mu1, mu2);
}

void BivariateNormalRand::Covariance(SquareMatrix<2> &matrix) const
{
    matrix(0, 0) = sigma1 * sigma1;
    matrix(0, 1) = matrix(1, 0) = rho * sigma1 * sigma2;
    matrix(1, 1) = sigma2 * sigma2;
}

double BivariateNormalRand::Correlation() const
{
    return rho;
}

void BivariateNormalRand::getFirstMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const
{
    distribution = X;
}

void BivariateNormalRand::getSecondMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const
{
    distribution = Y;
}

