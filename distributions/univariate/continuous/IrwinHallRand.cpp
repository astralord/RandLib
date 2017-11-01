#include "IrwinHallRand.h"

IrwinHallRand::IrwinHallRand(size_t number)
{
    SetNumber(number);
}

String IrwinHallRand::Name() const
{
    return "Irwin-Hall(" + toStringWithPrecision(GetNumber()) + ")";
}

void IrwinHallRand::SetNumber(int number)
{
    if (number <= 0)
        throw std::invalid_argument("Irwin-Hall distribution: number should be positive");
    n = number;
}

double IrwinHallRand::f(const double & x) const
{
    if (x < 0 || x > n)
        return 0.0;

    /// Check simplest cases
    if (n == 1)
        return 1.0;
    if (n == 2)
        return (x <= 1.0) ? x : 2.0 - x;
    /// General case
    double sum = 0.0;
    int last = std::floor(x);
    for (int i = 0; i <= last; ++i) {
        double add = (n - 1) * std::log(x - i);
        add -= RandMath::lfact(i);
        add -= RandMath::lfact(n - i);
        add = std::exp(add);
        sum += (i % 2) ? -add : add;
    }
    return n * sum;
}

double IrwinHallRand::logf(const double & x) const
{
    return std::log(f(x));
}

double IrwinHallRand::F(const double & x) const
{
    if (x <= 0)
        return 0.0;
    if (x >= n)
        return 1.0;

    /// Check simplest cases
    if (n == 1)
        return x;
    if (n == 2) {
        if (x <= 1.0)
            return 0.5 * x * x;
        double temp = 2.0 - x;
        return 1.0 - 0.5 * temp * temp;
    }
    /// General case
    double sum = 0.0;
    int last = std::floor(x);
    for (int i = 0; i <= last; ++i) {
        double add = n * std::log(x - i);
        add -= RandMath::lfact(i);
        add -= RandMath::lfact(n - i);
        add = std::exp(add);
        sum += (i % 2) ? -add : add;
    }
    return sum;
}

double IrwinHallRand::Variate() const
{
    double sum = 0.0;
    for (int i = 0; i != n; ++i)
        sum += U.Variate();
    return sum;
}

double IrwinHallRand::Mean() const
{
    return 0.5 * n;
}

double IrwinHallRand::Variance() const
{
    static constexpr double M_1_12 = 0.08333333333333;
    return n * M_1_12;
}

std::complex<double> IrwinHallRand::CFImpl(double t) const
{
    return std::pow(U.CF(t), n);
}

double IrwinHallRand::Median() const
{
    return 0.5 * n;
}

double IrwinHallRand::Mode() const
{
    return 0.5 * n;
}

double IrwinHallRand::Skewness() const
{
    return 0.0;
}

double IrwinHallRand::ExcessKurtosis() const
{
    return -1.2 / n;
}

