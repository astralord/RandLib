#pragma once

// Explicit instantiation declarations for the common types, active when RandLib is built with
// RANDLIB_PRECOMPILED=ON: the members are then compiled once, in src/precompiled.cpp, instead of
// in every translation unit that includes randlib.hpp.

#include <randlib/randlib.hpp>

namespace randlib
{

extern template class BivariateNormalRand<double>;
extern template class NormalInverseGammaRand<double>;
extern template class TrinomialRand<int>;
extern template class BetaPrimeRand<double>;
extern template class BetaDistribution<double>;
extern template class BetaRand<double>;
extern template class ArcsineRand<double>;
extern template class BaldingNicholsRand<double>;
extern template class CauchyRand<double>;
extern template class ContinuousDistribution<double>;
extern template class ExponentialRand<double>;
extern template class ExponentiallyModifiedGaussianRand<double>;
extern template class FisherFRand<double>;
extern template class FrechetRand<double>;
extern template class GeneralizedExtremeValueRand<double>;
extern template class GeneralizedParetoRand<double>;
extern template class GammaDistribution<double>;
extern template class FreeRateGammaDistribution<double>;
extern template class GammaRand<double>;
extern template class ChiSquaredRand<double>;
extern template class ErlangRand<double>;
extern template class GeneralGeometricStableDistribution<double>;
extern template class GeometricStableRand<double>;
extern template class GumbelRand<double>;
extern template class InverseGammaRand<double>;
extern template class InverseGaussianRand<double>;
extern template class IrwinHallRand<double>;
extern template class KolmogorovSmirnovRand<double>;
extern template class AsymmetricLaplaceDistribution<double>;
extern template class CenteredAsymmetricLaplaceRand<double>;
extern template class ShiftedAsymmetricLaplaceDistribution<double>;
extern template class AsymmetricLaplaceRand<double>;
extern template class LaplaceRand<double>;
extern template class LevyRand<double>;
extern template class LogNormalRand<double>;
extern template class LogisticRand<double>;
extern template class MarchenkoPasturRand<double>;
extern template class NakagamiDistribution<double>;
extern template class NakagamiRand<double>;
extern template class ChiRand<double>;
extern template class MaxwellBoltzmannRand<double>;
extern template class RayleighRand<double>;
extern template class NoncentralChiSquaredRand<double>;
extern template class NormalRand<double>;
extern template class ParetoRand<double>;
extern template class PlanckRand<double>;
extern template class RaisedCosineDistribution<double>;
extern template class RaisedCosineRand<double>;
extern template class RaabGreenRand<double>;
extern template class SechRand<double>;
extern template class StableDistribution<double>;
extern template class StableRand<double>;
extern template class HoltsmarkRand<double>;
extern template class LandauRand<double>;
extern template class StudentTRand<double>;
extern template class TriangularRand<double>;
extern template class UniformRand<double>;
extern template class WeibullRand<double>;
extern template class WignerSemicircleRand<double>;
extern template class CircularDistribution<double>;
extern template class VonMisesRand<double>;
extern template class WrappedExponentialRand<double>;
extern template class BetaBinomialRand<int>;
extern template class BinomialDistribution<int>;
extern template class BinomialRand<int>;
extern template class CategoricalRand<int>;
extern template class DiscreteDistribution<int>;
extern template class GeometricRand<int>;
extern template class HyperGeometricRand<int>;
extern template class LogarithmicRand<int>;
extern template class NegativeHyperGeometricRand<int>;
extern template class PoissonRand<int>;
extern template class SkellamRand<int>;
extern template class UniformDiscreteRand<int>;
extern template class YuleRand<int>;
extern template class ZetaRand<int>;
extern template class ZipfRand<int>;
extern template class NegativeBinomialRand<int, int>;
extern template class NegativeBinomialRand<int, double>;

} // namespace randlib
