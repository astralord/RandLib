#include "DiscreteBivariateDistribution.h"
#include "../univariate/discrete/BinomialRand.h"

template class DiscreteBivariateDistribution< BinomialRand<int>, BinomialRand<int>, int >;
template class DiscreteBivariateDistribution< BinomialRand<long int>, BinomialRand<long int>, long int >;
template class DiscreteBivariateDistribution< BinomialRand<long long int>, BinomialRand<long long int>, long long int >;
