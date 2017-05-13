#include "KumaraswamyRand.h"
#include "ParetoRand.h"
#include "UniformRand.h"

KumaraswamyRand::KumaraswamyRand(double shape1, double shape2)
{
    SetShapes(shape1, shape2);
}

std::string KumaraswamyRand::Name() const
{
    return "Kumaraswamy(" + toStringWithPrecision(GetFirstShape()) + ", "
                          + toStringWithPrecision(GetSecondShape()) + ")";
}

void KumaraswamyRand::SetShapes(double shape1, double shape2)
{
    a = (shape1 > 0) ? shape1 : 1.0;
    logA = std::log(a);
    b = (shape2 > 0) ? shape2 : 1.0;
    logB = std::log(b);
}

double KumaraswamyRand::f(const double & x) const
{
    if (x < 0.0 || x > 1.0)
        return 0.0;

    /// Deal with boundary cases
    if (x == 0) {
        if (a == 1)
            return a * b;
        return (a > 1) ? 0 : INFINITY;
    }
    if (x == 1) {
        if (b == 1)
            return a * b;
        return (b > 1) ? 0 : INFINITY;
    }
    return std::exp(logf(x));
}

double KumaraswamyRand::logf(const double & x) const
{
    if (x < 0.0 || x > 1.0)
        return -INFINITY;

    /// Deal with boundary cases
    if (x == 0) {
        if (a == 1)
            return logA + logB;
        return (a > 1) ? -INFINITY : INFINITY;
    }
    if (x == 1) {
        if (b == 1)
            return logA + logB;
        return (b > 1) ? -INFINITY : INFINITY;
    }

    double logX = std::log(x);
    double y = RandMath::log1mexp(a * logX);
    y *= b - 1;
    y += (a - 1) * logX;
    return logA + logB + y;
}

double KumaraswamyRand::F(const double & x) const
{
    if (x <= 0.0)
        return 0.0;
    if (x >= 1.0)
        return 1.0;
    double y = a * std::log(x);
    y = RandMath::log1mexp(y);
    return -std::expm1(b * y);
}

double KumaraswamyRand::S(const double & x) const
{
    if (x <= 0.0)
        return 1.0;
    if (x >= 1.0)
        return 0.0;
    double logX = std::log(x);
    double y = -std::expm1(a * logX); /// 1 - x^a
    return std::pow(y, b);
}

double KumaraswamyRand::Variate() const
{
    double X = 1.0 / ParetoRand::Variate(b, 1);
    return (X < 1e-5) ? std::exp(std::log1p(-X) / a) : std::pow(1.0 - X, 1.0 / a);
}

void KumaraswamyRand::Sample(std::vector<double> &outputData) const
{
    ParetoRand X(b);
    X.Sample(outputData);
    for (double & var : outputData) {
        var = 1.0 / var;
        if (var < 1e-5)
            var = (var < 1e-5) ? std::exp(std::log1p(-var) / a) : std::pow(1.0 - var, 1.0 / a);
    }
}

double KumaraswamyRand::Mean() const
{
    return Moment(1);
}

double KumaraswamyRand::Variance() const
{
    double m1 = Moment(1), m2 = Moment(2);
    return m2 - m1 * m1;
}

double KumaraswamyRand::Median() const
{
    return std::pow(-std::expm1(-M_LN2 / b), 1.0 / a);
}

double KumaraswamyRand::Mode() const
{
    if (a > 1)
        return (b > 1) ? std::pow((a - 1) / (a * b - 1), 1.0 / a) : 1.0;
    return (b > 1) ? 0.0 : (a > b);
}

double KumaraswamyRand::Skewness() const
{
    double m1 = Moment(1), m2 = Moment(2), m3 = Moment(3);
    double denominator = std::pow(m2 - m1 * m1, 1.5);
    double skewness = m3 - 3 * m1 * m2 + 2 * std::pow(m1, 3);
    return skewness / denominator;
}

double KumaraswamyRand::ExcessKurtosis() const
{
    double m1 = Moment(1), m2 = Moment(2), m3 = Moment(3), m4 = Moment(4);
    double var = m2 - m1 * m1;
    double denominator = var * var;
    double kurt = m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 - 3 * std::pow(m1, 4);
    return kurt / denominator;
}

double KumaraswamyRand::Moment(int n) const
{
    return b * RandMath::beta((a + n) / a, b);
}

double KumaraswamyRand::quantileImpl(double p) const
{
    if (p < 1e-5) {
        double x = std::log1p(-p);
        x = RandMath::log1mexp(x / b);
        return std::exp(x / a);
    }
    double x = std::pow(1.0 - p, 1.0 / b);
    return std::pow(1.0 - x, 1.0 / a);
}

double KumaraswamyRand::quantileImpl1m(double p) const
{
    double x = std::pow(p, 1.0 / b);
    if (x < 1e-5) {
        x = std::log1p(-x);
        return std::exp(x / a);
    }
    return std::pow(1.0 - x, 1.0 / a);
}
