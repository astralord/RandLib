#include "ProbabilityDistribution.h"
#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision

double2d & double2d::operator=(const double2d & other)
{
    x = other.x;
    y = other.y;
    return *this;
}

template < typename T >
ProbabilityDistribution<T>::ProbabilityDistribution()
{
}

template < typename T >
std::string ProbabilityDistribution<T>::toStringWithPrecision(const double a_value, const int n)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

template < typename T >
void ProbabilityDistribution<T>::sample(QVector<T> &outputData) const
{
    for (T &var : outputData)
        var = variate();
}

template < typename T >
void ProbabilityDistribution<T>::cdf(const QVector<T> &x, QVector<double> &y) const
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

template class ProbabilityDistribution<double>;
template class ProbabilityDistribution<double2d>;
template class ProbabilityDistribution< std::vector<double> >;
