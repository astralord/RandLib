#include "IrwinHallRand.h"

template < typename RealType >
IrwinHallRand<RealType>::IrwinHallRand(size_t number)
{
    SetNumber(number);
}

template < typename RealType >
String IrwinHallRand<RealType>::Name() const
{
    return "Irwin-Hall(" + this->toStringWithPrecision(GetNumber()) + ")";
}

template < typename RealType >
void IrwinHallRand<RealType>::SetNumber(int number)
{
    if (number <= 0)
        throw std::invalid_argument("Irwin-Hall distribution: number should be positive");
    n = number;
}

template < typename RealType >
double IrwinHallRand<RealType>::f(const RealType & x) const
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

template < typename RealType >
double IrwinHallRand<RealType>::logf(const RealType & x) const
{
    return std::log(f(x));
}

template < typename RealType >
double IrwinHallRand<RealType>::F(const RealType & x) const
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

template < typename RealType >
RealType IrwinHallRand<RealType>::Variate() const
{
    RealType sum = 0.0;
    for (int i = 0; i != n; ++i)
        sum += U.Variate();
    return sum;
}

template < typename RealType >
void IrwinHallRand<RealType>::Reseed(unsigned long seed) const
{
    U.Reseed(seed);
}

template < typename RealType >
long double IrwinHallRand<RealType>::Mean() const
{
    return 0.5 * n;
}

template < typename RealType >
long double IrwinHallRand<RealType>::Variance() const
{
    static constexpr long double M_1_12 = 0.08333333333333l;
    return n * M_1_12;
}

template < typename RealType >
std::complex<double> IrwinHallRand<RealType>::CFImpl(double t) const
{
    return std::pow(U.CF(t), n);
}

template < typename RealType >
RealType IrwinHallRand<RealType>::Median() const
{
    return 0.5 * n;
}

template < typename RealType >
RealType IrwinHallRand<RealType>::Mode() const
{
    return 0.5 * n;
}

template < typename RealType >
long double IrwinHallRand<RealType>::Skewness() const
{
    return 0.0l;
}

template < typename RealType >
long double IrwinHallRand<RealType>::ExcessKurtosis() const
{
    return -1.2 / n;
}

template class IrwinHallRand<float>;
template class IrwinHallRand<double>;
template class IrwinHallRand<long double>;
