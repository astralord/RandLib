#include "RaisedCosineRand.h"
#include "UniformRand.h"

RaisedCosineDistribution::RaisedCosineDistribution(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

void RaisedCosineDistribution::SetLocation(double location)
{
    mu = location;
}

void RaisedCosineDistribution::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Raised-Cosine distribution: scale should be positive");
    s = scale;
    s_pi = s * M_1_PI;
    log2S = std::log(2 * s);
}

double RaisedCosineDistribution::f(const double & x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -M_PI || xAdj >= M_PI)
        return 0.0;
    double y = std::cos(xAdj) + 1.0;
    return 0.5 * y / s;
}

double RaisedCosineDistribution::logf(const double & x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -M_PI || xAdj >= M_PI)
        return -INFINITY;
    double y = std::cos(xAdj);
    y = std::log1p(y);
    return y - log2S;
}

double RaisedCosineDistribution::F(const double & x) const
{
    double xAdj = (x - mu) / s;
    if (xAdj <= -1)
        return 0.0;
    if (xAdj >= 1)
        return 1.0;
    double y = std::sin(xAdj * M_PI);
    y *= M_1_PI;
    y += xAdj + 1;
    return 0.5 * y;
}

double RaisedCosineDistribution::S(const double & x) const
{
    double xAdj = (x - mu) / s;
    if (xAdj <= -1)
        return 1.0;
    if (xAdj >= 1)
        return 0.0;
    double y = std::sin(xAdj * M_PI);
    y *= M_1_PI;
    y += xAdj;
    return 0.5 - 0.5 * y;
}

double RaisedCosineDistribution::StandardVariate()
{
    /// p. 160. Non-Uniform Random Variate Generation. Luc Devroye
    double X = UniformRand::Variate(-M_PI_2, M_PI_2);
    double XSq = X * X;
    double U = UniformRand::StandardVariate();
    U += U;
    int a = 0, b = -1;
    double W = 0.0, V = 1.0;
    int iter = 0;
    do {
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W += V;
        if (U >= W)
            return X;
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W -= V;
        if (U <= W)
        {
            if (X == 0.0)
                return 0.0;
            return (X > 0.0) ? M_PI - X : -M_PI - X;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN;
}

double RaisedCosineDistribution::Variate() const
{
    return mu + s_pi * StandardVariate();
}

double RaisedCosineDistribution::Mean() const
{
    return mu;
}

double RaisedCosineDistribution::Variance() const
{
    static constexpr double coef = 1.0 / 3 - 2.0 / M_PI_SQ;
    return s * s * coef;
}

std::complex<double> RaisedCosineDistribution::CFImpl(double t) const
{
    double st = s * t;
    double numerator = M_PI_SQ * std::sin(st);
    double denominator = st * (M_PI_SQ - st * st);
    std::complex<double> y(0.0, mu * t);
    y = std::exp(y);
    return numerator / denominator * y;
}

double RaisedCosineDistribution::Median() const
{
    return mu;
}

double RaisedCosineDistribution::Mode() const
{
    return mu;
}

double RaisedCosineDistribution::Skewness() const
{
    return 0.0;
}

double RaisedCosineDistribution::ExcessKurtosis() const
{
    static constexpr double numerator = 1.2 * (90.0 - M_PI_SQ * M_PI_SQ);
    static constexpr double denominator = M_PI_SQ - 6.0;
    static constexpr double y = numerator / (denominator * denominator);
    return y;
}

String RaisedCosineRand::Name() const
{
    return "Raised cosine(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

String RaabGreenRand::Name() const
{
    return "Raab Green";
}
