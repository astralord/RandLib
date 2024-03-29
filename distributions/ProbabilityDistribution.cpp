#include "ProbabilityDistribution.h"
#include <sstream>
#include <iomanip>

template < typename T >
thread_local RandGenerator ProbabilityDistribution<T>::staticRandGenerator;

template < typename T >
ProbabilityDistribution<T>::ProbabilityDistribution()
{
}

template < typename T >
String ProbabilityDistribution<T>::toStringWithPrecision(const double a_value, const int n) const
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

template < typename T >
void ProbabilityDistribution<T>::CumulativeDistributionFunction(const std::vector<T> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = this->F(x[i]);
}

template < typename T >
double ProbabilityDistribution<T>::S(const T &x) const
{
    return 1.0 - this->F(x);
}

template < typename T >
void ProbabilityDistribution<T>::SurvivalFunction(const std::vector<T> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = this->S(x[i]);
}

template < typename T >
void ProbabilityDistribution<T>::Sample(std::vector<T> &outputData) const
{
    for (T &var : outputData)
        var = this->Variate();
}

template < typename T >
void ProbabilityDistribution<T>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
}

template < typename T >
constexpr char ProbabilityDistribution<T>::POSITIVITY_VIOLATION[];
template < typename T >
constexpr char ProbabilityDistribution<T>::NON_NEGATIVITY_VIOLATION[];
template < typename T >
constexpr char ProbabilityDistribution<T>::UPPER_LIMIT_VIOLATION[];
template < typename T >
constexpr char ProbabilityDistribution<T>::LOWER_LIMIT_VIOLATION[];

template < typename T >
String ProbabilityDistribution<T>::fitErrorDescription(ProbabilityDistribution::FIT_ERROR_TYPE fet, const String &explanation)
{
    String error = this->Name() + ": ";
    switch (fet) {
    case WRONG_SAMPLE:
        error += "Sample couldn't be generated by this distribution. ";
        break;
    case NOT_APPLICABLE:
        error += "Method cannot be applied here. ";
        break;
    case WRONG_RETURN:
        error += "Method returns invalid parameters. ";
        break;
    case TOO_FEW_ELEMENTS:
        error += "Sample is too small. ";
        break;
    case WRONG_LEVEL:
        error += "Significance level should be positive and smaller than 1. ";
        break;
    case UNDEFINED_ERROR:
    default:
        error += "Unknown type of error. ";
    }
    return error + explanation;
}

/// Univariate discrete
template class ProbabilityDistribution<int>;
template class ProbabilityDistribution<long int>;
template class ProbabilityDistribution<long long int>;

/// Univariate continuous
template class ProbabilityDistribution<float>;
template class ProbabilityDistribution<double>;
template class ProbabilityDistribution<long double>;

/// Bivariate discrete
template class ProbabilityDistribution< Pair<int> >;
template class ProbabilityDistribution< Pair<long int> >;
template class ProbabilityDistribution< Pair<long long int> >;

/// Bivariate continuous
template class ProbabilityDistribution< Pair<float> >;
template class ProbabilityDistribution< Pair<double> >;
template class ProbabilityDistribution< Pair<long double> >;
