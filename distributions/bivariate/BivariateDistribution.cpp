#include "BivariateDistribution.h"
#include "../univariate/discrete/BinomialRand.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"
#include "../univariate/continuous/InverseGammaRand.h"

template < class T1, class T2, typename T >
void BivariateDistribution<T1, T2, T>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    X.Reseed(seed + 1);
    Y.Reseed(seed + 2);
}

template < class T1, class T2, typename T >
LongDoublePair BivariateDistribution<T1, T2, T>::Mean() const
{
    return std::make_pair(X.Mean(), Y.Mean());
}

template < class T1, class T2, typename T >
LongDoubleTriplet BivariateDistribution<T1, T2, T>::Covariance() const
{
    double var1 = X.Variance();
    double var2 = Y.Variance();
    double corr = Correlation() * var1 * var2;
    return std::make_tuple(var1, corr, var2);
}

template < class T1, class T2, typename T >
std::pair<T1, T2> BivariateDistribution<T1, T2, T>::GetMarginalDistributions() const
{
    return std::make_pair(X, Y);
}

template class BivariateDistribution<NormalRand<float>, NormalRand<float>, float>;
template class BivariateDistribution<NormalRand<double>, NormalRand<double>, double>;
template class BivariateDistribution<NormalRand<long double>, NormalRand<long double>, long double>;

template class BivariateDistribution<StudentTRand<float>, InverseGammaRand<float>, float>;
template class BivariateDistribution<StudentTRand<double>, InverseGammaRand<double>, double>;
template class BivariateDistribution<StudentTRand<long double>, InverseGammaRand<long double>, long double>;

template class BivariateDistribution<BinomialRand<int>, BinomialRand<int>, int>;
template class BivariateDistribution<BinomialRand<long int>, BinomialRand<long int>, long int>;
template class BivariateDistribution<BinomialRand<long long int>, BinomialRand<long long int>, long long int>;
