#include "ContinuousBivariateDistribution.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"
#include "../univariate/continuous/InverseGammaRand.h"

template class ContinuousBivariateDistribution<NormalRand<float>, NormalRand<float>, float>;
template class ContinuousBivariateDistribution<NormalRand<double>, NormalRand<double>, double>;
template class ContinuousBivariateDistribution<NormalRand<long double>, NormalRand<long double>, long double>;

template class ContinuousBivariateDistribution<StudentTRand<float>, InverseGammaRand<float>, float>;
template class ContinuousBivariateDistribution<StudentTRand<double>, InverseGammaRand<double>, double>;
template class ContinuousBivariateDistribution<StudentTRand<long double>, InverseGammaRand<long double>, long double>;
