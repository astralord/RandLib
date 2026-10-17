#pragma once

#include <tuple>
#include <utility>

#if defined(RANDLIB_DISABLE_DEPRECATED_FIT_WARNINGS)
#define RANDLIB_DEPRECATED_FIT_MEMBER
#else
#define RANDLIB_DEPRECATED_FIT_MEMBER                                                                                  \
    [[deprecated("Use the non-mutating randlib::fit, randlib::fit_moments, or randlib::fit_bayes API")]]
#endif

#if defined(__clang__)
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN                                                                       \
    _Pragma("clang diagnostic push") _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"")
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END _Pragma("clang diagnostic pop")
#elif defined(__GNUC__)
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN                                                                       \
    _Pragma("GCC diagnostic push") _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END _Pragma("GCC diagnostic pop")
#elif defined(_MSC_VER)
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN __pragma(warning(push)) __pragma(warning(disable : 4996))
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END __pragma(warning(pop))
#else
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN
#define RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END
#endif

namespace randlib
{

template <typename T>
using Pair = std::pair<T, T>;

template <typename T>
using Triplet = std::tuple<T, T, T>;

using IntPair = Pair<int>;
using DoublePair = Pair<double>;
using LongDoublePair = Pair<long double>;

using IntTriplet = Triplet<int>;
using DoubleTriplet = Triplet<double>;
using LongDoubleTriplet = Triplet<long double>;

class not_implemented;
template <typename T>
class ProbabilityDistribution;
template <class T1, class T2, typename T>
class BivariateDistribution;
template <typename RealType = double>
class BivariateNormalRand;
template <class T1, class T2, typename RealType>
class ContinuousBivariateDistribution;
template <class T1, class T2, typename IntType>
class DiscreteBivariateDistribution;
template <typename RealType = double>
class NormalInverseGammaRand;
template <typename IntType = int>
class TrinomialRand;
template <typename T, typename P>
class ExponentialFamily;
template <typename T>
class UnivariateDistribution;
template <typename RealType = double>
class BetaPrimeRand;
template <typename RealType = double>
class BetaDistribution;
template <typename RealType = double>
class BetaRand;
template <typename RealType = double>
class ArcsineRand;
template <typename RealType = double>
class BaldingNicholsRand;
template <typename RealType = double>
class CauchyRand;
template <typename RealType>
class ContinuousDistribution;
class DegenerateRand;
template <typename RealType = double>
class ExponentialRand;
template <typename RealType = double>
class ExponentiallyModifiedGaussianRand;
template <typename RealType = double>
class FisherFRand;
template <typename RealType = double>
class FrechetRand;
template <typename RealType = double>
class GeneralizedExtremeValueRand;
template <typename RealType = double>
class GeneralizedParetoRand;
template <typename RealType = double>
class GammaDistribution;
template <typename RealType = double>
class FreeRateGammaDistribution;
template <typename RealType = double>
class GammaRand;
template <typename RealType = double>
class ChiSquaredRand;
template <typename RealType = double>
class ErlangRand;
template <typename RealType = double>
class GeneralGeometricStableDistribution;
template <typename RealType = double>
class GeometricStableRand;
template <typename RealType = double>
class GumbelRand;
template <typename RealType = double>
class InverseGammaRand;
template <typename RealType = double>
class InverseGaussianRand;
template <typename RealType = double>
class IrwinHallRand;
template <typename RealType = double>
class KolmogorovSmirnovRand;
template <typename RealType = double>
class AsymmetricLaplaceDistribution;
template <typename RealType = double>
class CenteredAsymmetricLaplaceRand;
template <typename RealType = double>
class ShiftedAsymmetricLaplaceDistribution;
template <typename RealType = double>
class AsymmetricLaplaceRand;
template <typename RealType = double>
class LaplaceRand;
template <typename RealType = double>
class LevyRand;
template <typename RealType = double>
class LogNormalRand;
template <typename RealType = double>
class LogisticRand;
template <typename RealType = double>
class LogStableRand;
template <typename RealType = double>
class MarchenkoPasturRand;
template <typename RealType = double>
class NakagamiDistribution;
template <typename RealType = double>
class NakagamiRand;
template <typename RealType = double>
class ChiRand;
template <typename RealType = double>
class MaxwellBoltzmannRand;
template <typename RealType = double>
class RayleighRand;
template <typename RealType = double>
class NoncentralChiSquaredRand;
template <typename RealType = double>
class NoncentralTRand;
template <typename RealType = double>
class NormalRand;
template <typename RealType = double>
class ParetoRand;
template <typename RealType = double>
class PlanckRand;
template <typename RealType = double>
class RaisedCosineDistribution;
template <typename RealType = double>
class RaisedCosineRand;
template <typename RealType = double>
class RaabGreenRand;
template <typename RealType = double>
class RiceRand;
template <typename RealType = double>
class SechRand;
template <typename RealType = double>
class StableDistribution;
template <typename RealType = double>
class StableRand;
template <typename RealType = double>
class StableCountRand;
template <typename RealType = double>
class HoltsmarkRand;
template <typename RealType = double>
class LandauRand;
template <typename RealType = double>
class StudentTRand;
template <typename RealType = double>
class TriangularRand;
template <typename RealType = double>
class UniformRand;
template <typename RealType = double>
class WeibullRand;
template <typename RealType = double>
class WignerSemicircleRand;
template <typename RealType = double>
class CircularDistribution;
template <typename RealType = double>
class VonMisesRand;
template <typename RealType = double>
class WrappedExponentialRand;
class BernoulliRand;
template <typename IntType = int>
class BetaBinomialRand;
template <typename IntType = int>
class BinomialDistribution;
template <typename IntType = int>
class BinomialRand;
template <typename IntType = int>
class CategoricalRand;
template <typename IntType>
class DiscreteDistribution;
template <typename IntType = int>
class GeometricRand;
template <typename IntType = int>
class HyperGeometricRand;
template <typename IntType = int>
class LogarithmicRand;
template <typename IntType = int, typename T = double>
class NegativeBinomialDistribution;
template <typename IntType = int, typename T = double>
class NegativeBinomialRand;
template <typename IntType = int>
class NegativeHyperGeometricRand;
template <typename IntType = int>
class PoissonRand;
class RademacherRand;
template <typename IntType = int>
class SkellamRand;
template <typename IntType = int>
class UniformDiscreteRand;
template <typename IntType = int>
class YuleRand;
template <typename IntType = int>
class ZetaRand;
template <typename IntType = int>
class ZipfRand;
class CantorRand;
class SingularDistribution;

} // namespace randlib
