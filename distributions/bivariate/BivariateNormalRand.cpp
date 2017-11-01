#include "BivariateNormalRand.h"


BivariateNormalRand::BivariateNormalRand(double location1, double location2, double scale1, double scale2, double correlation)
{
    SetLocations(location1, location2);
    SetCovariance(scale1, scale2, correlation);
}

String BivariateNormalRand::Name() const
{
    return "Bivariate Normal( (" + toStringWithPrecision(GetFirstLocation()) + ", "
                               + toStringWithPrecision(GetSecondLocation()) + "), ("
                               + toStringWithPrecision(GetFirstScale()) + ", "
                               + toStringWithPrecision(GetSecondScale()) + ", "
                               + toStringWithPrecision(GetCorrelation()) + ") )";
}

void BivariateNormalRand::SetLocations(double location1, double location2)
{
    mu1 = location1;
    mu2 = location2;
    X.SetLocation(mu1);
    Y.SetLocation(mu2);
}

void BivariateNormalRand::SetCovariance(double scale1, double scale2, double correlation)
{
    if (scale1 <= 0.0 || scale2 <= 0.0)
        throw std::invalid_argument("Scales of Bivariate-Normal distribution should be positive");
    if (std::fabs(correlation) > 1.0)
        throw std::invalid_argument("Correlation of Bivariate-Normal distribution should be in interval [-1, 1]");
    sigma1 = scale1;
    sigma2 = scale2;
    rho = correlation;

    X.SetScale(sigma1);
    Y.SetScale(sigma2);

    sqrt1mroSq = 0.5 * std::log1p(-rho * rho);
    pdfCoef = X.GetLogScale() + Y.GetLogScale() + sqrt1mroSq + M_LN2 + M_LNPI;
    sqrt1mroSq = std::exp(sqrt1mroSq);
}

double BivariateNormalRand::f(const DoublePair &point) const
{
    double xAdj = (point.first - mu1) / sigma1, yAdj = (point.second - mu2) / sigma2;
    if (rho == 1.0) /// f(x, y) = δ(xAdj - yAdj)
        return (xAdj - yAdj == 0) ? INFINITY : 0.0;
    if (rho == -1.0) /// f(x, y) = δ(xAdj + yAdj)
        return (xAdj + yAdj == 0) ? INFINITY : 0.0;
    return std::exp(logf(point));
}

double BivariateNormalRand::logf(const DoublePair &point) const
{
    if (rho == 0.0) /// log(f(x, y)) = log(f(x)) + log(f(y))
        return X.logf(point.first) + Y.logf(point.second);
    double xAdj = (point.first - mu1) / sigma1, yAdj = (point.second - mu2) / sigma2;
    if (rho == 1.0) /// log(f(x, y)) = log(δ(xAdj - yAdj))
        return (xAdj - yAdj == 0) ? INFINITY : -INFINITY;
    if (rho == -1.0) /// log(f(x, y)) = log(δ(xAdj + yAdj))
        return (xAdj + yAdj == 0) ? INFINITY : -INFINITY;
    double z = xAdj * xAdj - 2 * rho * xAdj * yAdj + yAdj * yAdj;
    return -(pdfCoef + 0.5 * z / (sqrt1mroSq * sqrt1mroSq));
}

double BivariateNormalRand::F(const DoublePair &point) const
{
    if (rho == 0.0) /// F(x, y) = F(x)F(y)
        return X.F(point.first) * Y.F(point.second);
    /// Unnumbered equation between (3) and (4) in Section 2.2 of Genz (2004),
    /// integrating in terms of theta between asin(ρ) and +/- π/2
    double p1, p2 = 0.0;
    double x = point.first, y = point.second;
    double xAdj = (x - mu1) / sigma1, yAdj = (y - mu2) / sigma2;
    if (rho > 0)
        p1 = (xAdj < yAdj) ? X.F(x) : Y.F(y);
    else
        p1 = std::max(X.F(x) - Y.F(-y), 0.0);
    if (std::fabs(rho) < 1.0) {
        double lowLimit = std::asin(rho);
        double highLimit = RandMath::sign(rho) * M_PI_2;
        p2 = RandMath::integral([this, xAdj, yAdj] (double theta) {
            /// Integrand is exp(-(x^2 + y^2 - 2xysin(θ)) / (2cos(θ)^2))
            double cosTheta = std::cos(theta);
            double tanTheta = std::tan(theta);
            double integrand = xAdj * tanTheta - yAdj / cosTheta;
            integrand *= integrand;
            integrand += xAdj * xAdj;
            integrand = std::exp(-0.5 * integrand);
            return integrand;
        }, lowLimit, highLimit);
    }
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

double BivariateNormalRand::Correlation() const
{
    return rho;
}
