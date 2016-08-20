#include "WaldRand.h"
#include "NormalRand.h"
#include "UniformRand.h"

WaldRand::WaldRand(double mean, double shape)
{
    setParameters(mean, shape);
}

std::string WaldRand::name() const
{
    return "Wald(" + toStringWithPrecision(getMean()) + ", " + toStringWithPrecision(getShape()) + ")";
}

void WaldRand::setParameters(double mean, double shape)
{
    mu = mean;
    if (mu <= 0)
        mu = 1.0;
        
    lambda = shape;
    if (lambda <= 0)
        lambda = 1.0;
    pdfCoef = 0.5 * std::log(0.5 * lambda * M_1_PI);
    cdfCoef = std::exp(2 * lambda / mu);
}

double WaldRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double xInv = 1.0 / x;
    double y = xInv * std::sqrt(xInv);
    double z = (x - mu);
    z *= z;
    z *= -.5 * lambda * xInv / (mu * mu);
    z = std::exp(z + pdfCoef);
    return y * z;
}

double WaldRand::F(double x) const
{
    if (x <= 0)
        return 0;
    double a = std::sqrt(0.5 * lambda / x);
    double b = a * x / mu;
    double y = 1 + std::erf(b - a);
    double z = std::erfc(b + a);
    return 0.5 * (y + cdfCoef * z);
}

double WaldRand::variate() const
{
    double y = NormalRand::standardVariate();
    y *= y;
    double my = mu * y;
    double x = 4 * lambda + my;
    x *= my;
    x = std::sqrt(x);
    x = my - x;
    x *= .5 / lambda;
    ++x;
    return (UniformRand::standardVariate() <= 1.0 / (1 + x)) ? mu * x : mu / x;
}

std::complex<double> WaldRand::CF(double t) const
{
    if (t == 0)
        return 1;
    double im = mu * mu;
    im *= t / lambda;
    std::complex<double> y(1, -im - im);
    y = 1.0 - std::sqrt(y);
    y *= lambda / mu;
    return std::exp(y);
}

double WaldRand::Mode() const
{
    double aux = 1.5 * mu / lambda;
    double mode = 1 + aux * aux;
    mode = std::sqrt(mode);
    mode -= aux;
    return mu * mode;
}

double WaldRand::Skewness() const
{
    return 3 * std::sqrt(mu / lambda);
}

double WaldRand::ExcessKurtosis() const
{
    return 15 * mu / lambda;
}
