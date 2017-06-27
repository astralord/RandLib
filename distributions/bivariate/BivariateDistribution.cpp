#include "BivariateDistribution.h"
#include "../univariate/discrete/BinomialRand.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"
#include "../univariate/continuous/InverseGammaRand.h"
#include "../univariate/discrete/BinomialRand.h"

template < class T1, class T2, typename T >
DoublePair BivariateDistribution<T1, T2, T>::Mean() const
{
    return std::make_pair(X.Mean(), Y.Mean());
}

template < class T1, class T2, typename T >
DoubleTriplet BivariateDistribution<T1, T2, T>::Covariance() const
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

template class BivariateDistribution<NormalRand, NormalRand, DoublePair>;
template class BivariateDistribution<StudentTRand, InverseGammaRand, DoublePair>;
template class BivariateDistribution<BinomialRand, BinomialRand, IntPair>;
