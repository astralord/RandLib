#ifndef RANDOMGENERATORS_H
#define RANDOMGENERATORS_H

#include <RandomVariable.h>

/// DISCRETE
#include <discrete/UniformDiscreteRand.h>
#include <discrete/PoissonRand.h>

/// CONTINUOUS
/// Most used
#include <continuous/UniformRand.h>
#include <continuous/NormalRand.h>
#include <continuous/ExponentialRand.h>
#include <continuous/GammaRand.h>
/// Fairly common
#include <continuous/CauchyRand.h>
#include <continuous/LevyRand.h>
#include <continuous/LaplaceRand.h>
#include <continuous/LogNormalRand.h>
#include <continuous/ChiSquaredRand.h>
#include <continuous/ErlangRand.h>
#include <continuous/BetaRand.h>
#include <continuous/StudentTRand.h>
#include <continuous/FisherSnedecorRand.h>
#include <continuous/ParetoRand.h>
#include <continuous/RayleighRand.h>
/// Advanced
#include <continuous/StableRand.h>


#endif // RANDOMGENERATORS_H
