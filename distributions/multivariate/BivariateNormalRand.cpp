#include "BivariateNormalRand.h"


BivariateNormalRand::BivariateNormalRand(const DoublePair &location, const DoubleTriplet &covariance)
{
    SetLocation(location);
    SetScale(covariance);
}

BivariateNormalRand::BivariateNormalRand(double location1, double location2, double scale1, double correlation, double scale2)
{
    SetLocation(location1, location2);
    SetScale(scale1, correlation, scale2);
}

std::string BivariateNormalRand::Name() const
{
    return "Bivariate Normal( (" + toStringWithPrecision(GetFirstLocation()) + ", "
                               + toStringWithPrecision(GetSecondLocation()) + "), ("
                               + toStringWithPrecision(GetFirstScale()) + ", "
                               + toStringWithPrecision(GetSecondScale()) + ", "
                               + toStringWithPrecision(GetCorrelation()) + ") )";
}

void BivariateNormalRand::SetLocation(const DoublePair &location)
{
    mu1 = location.first;
    mu2 = location.second;
    X.SetLocation(mu1);
    Y.SetLocation(mu2);
}

void BivariateNormalRand::SetLocation(double location1, double location2)
{
    mu1 = location1;
    mu2 = location2;
    X.SetLocation(mu1);
    Y.SetLocation(mu2);
}

void BivariateNormalRand::SetScale(double scale1, double correlation, double scale2)
{
    sigma1 = scale1 > 0 ? scale1 : 1.0;
    sigma2 = scale2 > 0 ? scale2 : 1.0;
    rho = correlation;
    if (rho < 0 || rho > 1)
        rho = 0.0;

    X.SetScale(sigma1);
    Y.SetScale(sigma2);

    sqrt1mroSq = std::sqrt(1.0 - rho * rho);
    pdfCoef = 0.5 * M_1_PI / (sigma1 * sigma2 * sqrt1mroSq);
}

void BivariateNormalRand::SetScale(const DoubleTriplet &covariance)
{
    double var1, cov, var2;
    std::tie(var1, cov, var2) = covariance;
    double scale1 = std::sqrt(var1);
    double scale2 = std::sqrt(var2);
    double corr = cov * scale1 * scale2;
    SetScale(std::sqrt(var1), corr, std::sqrt(var2));
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

DoublePair BivariateNormalRand::Variate() const
{
    double Z1 = NormalRand::StandardVariate();
    double Z2 = NormalRand::StandardVariate();
    double x = mu1 + sigma1 * Z1;
    double y = mu2 + sigma2 * (rho * Z1 + sqrt1mroSq * Z2);
    return std::make_pair(x, y);
}

DoublePair BivariateNormalRand::Mean() const
{
    return std::make_pair(mu1, mu2);
}

DoubleTriplet BivariateNormalRand::Covariance() const
{
    double var1 = sigma1 * sigma1;
    double var2 = sigma2 * sigma2;
    double corr = rho * sigma1 * sigma2;
    return std::make_tuple(var1, corr, var2);
}

double BivariateNormalRand::Correlation() const
{
    return rho;
}

void BivariateNormalRand::GetMarginalDistributions(UnivariateProbabilityDistribution<double> &distribution1, UnivariateProbabilityDistribution<double> &distribution2) const
{
    distribution1 = X;
    distribution2 = Y;
}
