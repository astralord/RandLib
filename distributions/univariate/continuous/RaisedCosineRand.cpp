#include "RaisedCosineRand.h"
#include "UniformRand.h"

RaisedCosineRand::RaisedCosineRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string RaisedCosineRand::name()
{
    return "Raised cosine(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void RaisedCosineRand::setLocation(double location)
{
    mu = location;
}

void RaisedCosineRand::setScale(double scale)
{
    s = scale;
    if (s <= 0)
        s = M_PI;
    s_pi = s * M_1_PI;
}

double RaisedCosineRand::f(double x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -M_PI || xAdj >= M_PI)
        return 0.0;
    double y = std::cos(xAdj);
    ++y;
    return y / (s + s);
}

double RaisedCosineRand::F(double x) const
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

double RaisedCosineRand::standardVariate()
{
    /// p. 160. Non-Uniform Random Variate Generation. Luc Devroye
    double X = UniformRand::variate(-M_PI_2, M_PI_2);
    double XSq = X * X;
    double U = UniformRand::standardVariate();
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
    } while (++iter < 1e9);
    return NAN;
}

double RaisedCosineRand::variate() const
{
    return mu + s_pi * standardVariate();
}

double RaisedCosineRand::Mean() const
{
    return mu;
}

double RaisedCosineRand::Variance() const
{
    double y = M_1_PI;
    y *= y;
    y += y;
    y += 1.0 / 3;
    return s * s * y;
}

std::complex<double> RaisedCosineRand::CF(double t) const
{
    double st = s * t;
    double numerator = M_PI_SQ * std::sin(st);
    double denominator = st * (M_PI_SQ - st * st);
    std::complex<double> y(0.0, mu * t);
    y = std::exp(y);
    return numerator / denominator * y;
}

double RaisedCosineRand::Median() const
{
    return mu;
}

double RaisedCosineRand::Mode() const
{
    return mu;
}

double RaisedCosineRand::Skewness() const
{
    return 0.0;
}

double RaisedCosineRand::ExcessKurtosis() const
{
    static constexpr double numerator = 1.2 * (90.0 - M_PI_SQ * M_PI_SQ);
    static constexpr double denominator = M_PI_SQ - 6.0;
    static constexpr double y = numerator / (denominator * denominator);
    return y;
}


std::string RaabGreenRand::name()
{
    return "Raab Green";
}
