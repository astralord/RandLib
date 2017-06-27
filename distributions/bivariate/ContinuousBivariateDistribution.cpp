#include "ContinuousBivariateDistribution.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"
#include "../univariate/continuous/InverseGammaRand.h"

template class ContinuousBivariateDistribution<NormalRand, NormalRand>;
template class ContinuousBivariateDistribution<StudentTRand, InverseGammaRand>;
