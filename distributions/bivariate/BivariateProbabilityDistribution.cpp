#include "BivariateProbabilityDistribution.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"
#include "../univariate/continuous/InverseGammaRand.h"

template < class T1, class T2 >
DoubleTriplet BivariateProbabilityDistribution<T1, T2>::Covariance() const
{
    double var1 = X.Variance();
    double var2 = Y.Variance();
    double corr = Correlation();
    return std::make_tuple(var1, corr, var2);
}

template < class T1, class T2 >
std::pair<T1, T2> BivariateProbabilityDistribution<T1, T2>::GetMarginalDistributions() const
{
    return std::make_pair(X, Y);
}

template class BivariateProbabilityDistribution<NormalRand, NormalRand>;
template class BivariateProbabilityDistribution<StudentTRand, InverseGammaRand>;
