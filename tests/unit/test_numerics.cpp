// Quadrature, root finding and optimization in randlib::math.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <numbers>
#include <span>
#include <vector>

using namespace randlib;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Numerical functions preserve long-double arithmetic", "[numerics]")
{
    using Real = long double;
    auto checkClose = [](Real actual, Real expected, Real relativeTolerance = 3e-16L) {
        CHECK(std::fabs(actual - expected) <= relativeTolerance * std::max(Real{1}, std::fabs(expected)));
    };

    STATIC_CHECK(std::same_as<decltype(math::softplus(Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::logBeta(Real{1}, Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::pgamma(Real{1}, Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::besselI(Real{0}, Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::factorial(Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::factorial(1)), long double>);
    STATIC_CHECK(std::same_as<decltype(math::doubleFactorial<Real>(1)), Real>);
    STATIC_CHECK(std::same_as<decltype(math::lfact<Real>(1)), Real>);

    checkClose(math::atan(Real{2}), std::atan(Real{2}));
    checkClose(math::softplus(Real{-2}), std::log1p(std::exp(Real{-2})));
    checkClose(math::log1mexp(Real{-2}), std::log1p(-std::exp(Real{-2})));
    checkClose(math::logexpm1l(Real{2}), std::log(std::expm1(Real{2})));
    checkClose(math::log2mexp(Real{-2}), std::log(Real{2} - std::exp(Real{-2})));

    constexpr Real inverseArgument = 0.375L;
    checkClose(math::erfinv(std::erf(inverseArgument)), inverseArgument, 2e-15L);
    checkClose(math::erfcinv(std::erfc(inverseArgument)), inverseArgument, 2e-15L);
    checkClose(math::erfcx(Real{2}), std::exp(Real{4}) * std::erfc(Real{2}), 2e-15L);
    checkClose(math::harmonicNumber(Real{1}, 4), Real{25} / 12, 1e-14L);

    checkClose(math::factorial(Real{10}), Real{3628800});
    checkClose(math::lfact<Real>(10), std::lgamma(Real{11}));
    checkClose(math::ldfact<Real>(9), std::log(Real{945}));
    checkClose(math::binom<Real>(20, 6), Real{38760}, 4e-15L);

    constexpr Real x = 2.25L;
    checkClose(math::digamma(x + 1), math::digamma(x) + 1 / x);
    checkClose(math::trigamma(x + 1), math::trigamma(x) - 1 / (x * x));
    checkClose(math::ibeta(Real{0.25}, Real{2}, Real{3}), Real{0.26171875}, 2e-15L);
    checkClose(math::pgamma(Real{1}, x), -std::expm1(-x), 2e-15L);
    checkClose(math::pgamma(Real{2}, x) + math::qgamma(Real{2}, x), Real{1}, 2e-15L);

    constexpr Real besselArgument = 1.25L;
    Real expectedBesselI = std::sqrt(Real{2} / (std::numbers::pi_v<Real> * besselArgument)) * std::sinh(besselArgument);
    checkClose(math::portable::besselI(Real{0.5}, besselArgument), expectedBesselI, 2e-15L);
    checkClose(math::portable::riemannZeta(Real{2}), std::numbers::pi_v<Real> * std::numbers::pi_v<Real> / 6, 2e-15L);

    Real w = math::W0Lambert(Real{3});
    checkClose(w * std::exp(w), Real{3}, 2e-15L);
    Real marcumP = math::MarcumP(Real{1}, Real{2}, Real{3});
    Real marcumQ = math::MarcumQ(Real{1}, Real{2}, Real{3});
    CHECK(marcumP >= 0);
    CHECK(marcumP <= 1);
    checkClose(marcumP + marcumQ, Real{1});

    constexpr Real tailArgument = 751.09573767168045L;
    CHECK_THAT(static_cast<double>(math::MarcumQ(Real{0.25}, Real{1.5}, tailArgument)), WithinRel(1e-300, 1e-12));
}

TEST_CASE("Adaptive quadrature on finite and infinite intervals", "[numerics]")
{
    CHECK_THAT(math::integrate([](double x) { return std::exp(x); }, 0, 1), WithinRel(std::numbers::e - 1, 1e-13));
    CHECK_THAT(math::integrate([](double x) { return 1 / std::sqrt(x); }, 0, 1, 1e-12, 4096), WithinRel(2.0, 1e-10));
    CHECK_THAT(math::integrate([](double x) { return std::exp(-x * x); }, -INFINITY, INFINITY),
               WithinRel(std::sqrt(std::numbers::pi), 1e-12));
    CHECK_THAT(math::integrate([](double x) { return std::exp(-x); }, 2, INFINITY), WithinRel(std::exp(-2.0), 1e-12));
    CHECK_THAT(math::integrate([](double x) { return 1 / (1 + x * x); }, -INFINITY, 0),
               WithinRel(std::numbers::pi / 2, 1e-12));
    CHECK(math::integrate([](double x) { return x; }, 1, 0) == -math::integrate([](double x) { return x; }, 0, 1));
    CHECK_THAT(static_cast<double>(math::integrate<long double>([](long double x) { return x * x; }, 0, 3)),
               WithinRel(9.0, 1e-15));
}

TEST_CASE("Adaptive quadrature stops when rounding errors dominate", "[numerics]")
{
    /// The integrand carries deterministic noise of relative size 1e-13, so no subdivision can reach 1e-15
    int evaluations = 0;
    const double value = math::integrate(
        [&evaluations](double x) {
            ++evaluations;
            const double noise = std::sin(1e7 * x) * 1e-13;
            return std::exp(x) * (1 + noise);
        },
        0, 1, 1e-15, 262144);
    CHECK_THAT(value, WithinRel(std::numbers::e - 1, 1e-12));
    CHECK(evaluations < 21 * 2 * 200);
}

TEST_CASE("Noncentral t series, densities and moments match high-precision values", "[numerics][noncentral-t]")
{
    const NoncentralTRand<double> distribution(30, 5);
    CHECK_THAT(distribution.F(-10), WithinRel(1.4196688792261714e-25, 1e-12));
    CHECK_THAT(distribution.f(-10), WithinRel(3.8127458982120124e-25, 1e-12));
    CHECK_THAT(distribution.F(-3), WithinRel(3.2626878724789141e-14, 1e-12));
    CHECK_THAT(distribution.f(-3), WithinRel(1.6499531574158244e-13, 1e-12));
    CHECK_THAT(distribution.f(1), WithinRel(0.0001552857352701685070297871, 1e-12));
    CHECK_THAT(distribution.Mean(), WithinRel(5.1294973500495874, 1e-14));
    CHECK_THAT(distribution.Variance(), WithinRel(1.5453997929771642, 1e-13));
    CHECK_THAT(distribution.Skewness(), WithinRel(0.43566262458486494, 1e-12));
    CHECK_THAT(distribution.ExcessKurtosis(), WithinRel(0.5001967695937044, 1e-11));
    CHECK(distribution.Quantile(distribution.F(0)) == 0);
    CHECK(distribution.Quantile1m(distribution.S(0)) == 0);
    const auto characteristicFunction = distribution.CF(0.5);
    CHECK_THAT(characteristicFunction.real(), WithinAbs(-0.68531584006586704, 2e-12));
    CHECK_THAT(characteristicFunction.imag(), WithinAbs(0.46250925896885466, 2e-12));
    const auto lowDegreeCf = NoncentralTRand<double>(1, 1).CF(0.5);
    CHECK_THAT(lowDegreeCf.real(), WithinAbs(0.47201952890831652, 2e-12));
    CHECK_THAT(lowDegreeCf.imag(), WithinAbs(0.31943044765159106, 2e-12));
    const auto negativeCf = NoncentralTRand<double>(5, -2).CF(0.5);
    CHECK_THAT(negativeCf.real(), WithinAbs(0.35042418454890825, 2e-12));
    CHECK_THAT(negativeCf.imag(), WithinAbs(-0.68761133161353897, 2e-12));

    const NoncentralTRand<double> reflected(30, -5);
    CHECK_THAT(distribution.F(-3), WithinRel(reflected.S(3), 1e-14));
    CHECK_THAT(distribution.f(-3), WithinRel(reflected.f(3), 1e-14));

    const NoncentralTRand<double> central(5, 0);
    const StudentTRand<double> student(5);
    for (double x : {-10.0, -1.0, 0.0, 2.0, 20.0}) {
        CHECK_THAT(central.f(x), WithinRel(student.f(x), 1e-14));
        CHECK_THAT(central.F(x), WithinAbs(student.F(x), 1e-14));
        CHECK_THAT(central.S(x), WithinAbs(student.S(x), 1e-14));
    }

    const NoncentralTRand<double> largeNoncentrality(30, 100);
    CHECK_THAT(largeNoncentrality.F(101), WithinRel(0.49626388164838351, 1e-13));
    CHECK_THAT(largeNoncentrality.S(101), WithinRel(0.50373611835161649, 1e-13));
    CHECK_THAT(largeNoncentrality.f(50), WithinRel(2.7944329641841315e-12, 1e-12));
    CHECK_THAT(largeNoncentrality.f(101), WithinRel(0.030247658482470235, 1e-12));
    CHECK_THAT(largeNoncentrality.S(2000), WithinRel(3.1437473004343311e-34, 1e-11));
    CHECK_THAT(largeNoncentrality.f(2000), WithinRel(4.7045030926584685e-36, 1e-11));
}

TEST_CASE("Tanh-sinh quadrature resolves endpoint singularities", "[numerics]")
{
    const auto betaHalf = [](double x, double endpointDistance) {
        const double leftDistance = endpointDistance < 0 ? -endpointDistance : x;
        const double rightDistance = endpointDistance > 0 ? endpointDistance : 1 - x;
        return 1 / std::sqrt(leftDistance * rightDistance);
    };
    CHECK_THAT(math::integrateTanhSinh(betaHalf, 0, 1), WithinRel(std::numbers::pi, 2e-13));
    CHECK_THAT(math::integrateTanhSinh([](double x) { return std::log(x) * std::log1p(-x); }, 0, 1),
               WithinRel(2 - std::numbers::pi * std::numbers::pi / 6, 2e-13));

    const auto endpointAware = [](double x, double endpointDistance) {
        if (x < 0.5)
            return std::sqrt(x / (1 - x * x));
        return std::sqrt(x / ((1 + x) * endpointDistance));
    };
    const double expected = 0.5 * std::tgamma(0.75) * std::tgamma(0.5) / std::tgamma(1.25);
    CHECK_THAT(math::integrateTanhSinh(endpointAware, 0, 1), WithinRel(expected, 2e-13));

    const long double result =
        math::integrateTanhSinh<long double>([](long double x) { return 1 / std::sqrt(x); }, 0, 1, 1e-15L);
    CHECK_THAT(static_cast<double>(result), WithinRel(2.0, 2e-15));
    CHECK(math::integrateTanhSinh([](double x) { return x; }, 1, 0) ==
          -math::integrateTanhSinh([](double x) { return x; }, 0, 1));
    CHECK_THROWS_AS(math::integrateTanhSinh([](double x) { return x; }, 0, INFINITY), std::invalid_argument);
}

TEST_CASE("Oscillatory quadrature accelerates semi-infinite Fourier integrals", "[numerics]")
{
    for (double frequency : {0.1, 1.0, 5.0}) {
        const auto actual = math::integrateOscillatory([](double x) { return std::exp(-x); }, frequency);
        const std::complex<double> expected = 1.0 / std::complex<double>(1.0, -frequency);
        CHECK_THAT(std::abs(actual - expected), WithinAbs(0.0, 1e-10));
    }

    const auto pareto = math::integrateOscillatory([](double x) { return 0.5 * std::pow(x, -1.5); }, 0.1, 1.0);
    CHECK_THAT(pareto.real(), WithinAbs(0.605333341794216409742774, 1e-8));
    CHECK_THAT(pareto.imag(), WithinAbs(0.2963660538362012702997561, 1e-8));
    CHECK_THROWS_AS(math::integrateOscillatory([](double x) { return x; }, 0.0), std::invalid_argument);

    /// A Gaussian bump many half-periods above the lower limit: empty half-periods do not end the walk
    /// when the tail mass is known
    const double center = 30, width = 0.1, frequency = 20;
    auto bump = [=](double x) {
        return std::exp(-0.5 * std::pow((x - center) / width, 2)) / (width * std::sqrt(2 * std::numbers::pi));
    };
    auto tail = [=](double x) { return 0.5 * std::erfc((x - center) / (width * std::numbers::sqrt2)); };
    const auto shifted = math::integrateOscillatory(bump, frequency, 0.0, tail);
    const std::complex<double> gaussian =
        std::polar(std::exp(-0.5 * std::pow(width * frequency, 2)), frequency * center);
    CHECK_THAT(std::abs(shifted - gaussian), WithinAbs(0.0, 1e-10));
}

TEST_CASE("Incomplete beta function keeps accuracy for rounded complements and large shapes", "[numerics]")
{
    /// x = 1 - 1e-30 rounds to 1; the logarithms carry the exact complement
    const double tiny = 1e-30;
    const double logBeta = math::logBeta(0.442, 0.0175);
    CHECK_THAT(math::ibeta(1.0 - tiny, 0.442, 0.0175, logBeta, std::log1p(-tiny), std::log(tiny)),
               WithinRel(0.71004226954528407701, 1e-14));
    /// Binomial(2e9, 1e-8) needs about 2e9 terms of the finite continued fractions
    const BinomialRand<int> binomial(2000000000, 1e-8);
    CHECK_THAT(binomial.F(17), WithinRel(0.29702839678536085955, 1e-11));
    CHECK_THAT(binomial.F(20), WithinRel(0.55909258423132520632, 1e-11));
    CHECK_THAT(binomial.F(20) + binomial.S(20), WithinAbs(1.0, 1e-14));
}

TEST_CASE("Marcum functions are complementary for large noncentrality", "[numerics]")
{
    /// Rice(50, 0.1) has halfNoncentrality 125000, beyond the small-argument series
    const RiceRand<double> rice(50, 0.1);
    CHECK_THAT(rice.F(50), WithinRel(0.49960105752012697831, 1e-8));
    CHECK_THAT(rice.S(50.8), WithinRel(6.2712828198160881433e-16, 1e-8));
    for (double x : {49.0, 49.95, 50.0, 50.3})
        CHECK_THAT(rice.F(x) + rice.S(x), WithinAbs(1.0, 4e-16));
    const NoncentralChiSquaredRand<double> chiSquared(10, 40000);
    CHECK_THAT(chiSquared.F(40010), WithinRel(0.5009972549393809418, 1e-9));
    CHECK_THAT(chiSquared.F(40010) + chiSquared.S(40010), WithinAbs(1.0, 4e-16));
}

TEST_CASE("Logarithms of Bessel functions beyond the range of double", "[numerics]")
{
    CHECK_THAT(math::logBesselI<double>(0, 1000), WithinRel(995.62730888986946467, 1e-15));
    CHECK_THAT(math::logBesselI<double>(268, 12.45), WithinRel(-743.89755522925462318, 1e-14));
    CHECK_THAT(math::logBesselK<double>(0, 800), WithinRel(-803.11667066365989547, 1e-15));
    CHECK_THAT(math::logBesselK<double>(-2.5, 3), WithinRel(-2.4762169313021237996, 1e-14));
}

TEST_CASE("Complex Lambert W function on the principal branch", "[numerics]")
{
    const auto unit = math::W0Lambert(std::complex<double>(0, -1));
    CHECK_THAT(std::abs(unit - std::complex<double>(0.37469902073711749361, -0.57641272303143528315)),
               WithinAbs(0.0, 1e-15));
    const auto large = math::W0Lambert(std::complex<double>(0, -1e15));
    CHECK_THAT(std::abs(large - std::complex<double>(31.100361073551926954, -1.5219002018854649007)),
               WithinAbs(0.0, 1e-13));
    for (std::complex<double> z : {std::complex<double>(0, 1e-12), {2, 3}, {-0.3, 0.1}, {1, 0}}) {
        const auto w = math::W0Lambert(z);
        CHECK_THAT(std::abs(w * std::exp(w) - z), WithinAbs(0.0, 1e-15 * std::abs(z)));
    }
}

TEST_CASE("Sample statistics merge chunks of ten thousand observations", "[numerics]")
{
    std::vector<double> twoBlocks(20000);
    for (std::size_t i = 0; i < twoBlocks.size(); ++i)
        twoBlocks[i] = (i < 10000) ? 0.0 : 1.0;
    const auto [mean, variance, skewness, excessKurtosis] = NormalRand<double>::GetSampleStatistics(twoBlocks);
    CHECK_THAT(static_cast<double>(mean), WithinRel(0.5, 1e-15));
    CHECK_THAT(static_cast<double>(variance), WithinRel(0.25, 1e-12));
    CHECK_THAT(static_cast<double>(skewness), WithinAbs(0.0, 1e-12));
    CHECK_THAT(static_cast<double>(excessKurtosis), WithinRel(-2.0, 1e-12));

    std::vector<double> ramp(25000);
    for (std::size_t i = 0; i < ramp.size(); ++i)
        ramp[i] = static_cast<double>(i);
    const double n = static_cast<double>(ramp.size());
    const auto [rampMean, rampVariance, rampSkewness, rampKurtosis] = NormalRand<double>::GetSampleStatistics(ramp);
    CHECK_THAT(static_cast<double>(rampMean), WithinRel((n - 1) / 2, 1e-15));
    CHECK_THAT(static_cast<double>(rampVariance), WithinRel((n * n - 1) / 12, 1e-12));
    CHECK_THAT(static_cast<double>(rampSkewness), WithinAbs(0.0, 1e-12));
    CHECK_THAT(static_cast<double>(rampKurtosis), WithinRel(-1.2 * (n * n + 1) / (n * n - 1), 1e-9));
}

TEST_CASE("TOMS 748 finds bracketed roots to full precision", "[numerics]")
{
    double root = math::findRoot([](double x) { return std::cos(x) - x; }, 0, 1);
    CHECK_THAT(root, WithinRel(0.73908513321516064166, 4e-16));
    double cubic = math::findRoot([](double x) { return x * x * x - 2; }, 0, 10);
    CHECK_THAT(cubic, WithinRel(std::cbrt(2.0), 4e-16));
    /// a steep function with a root near one end
    double steep = math::findRoot([](double x) { return std::exp(20 * x) - 1e-3; }, -5, 5);
    CHECK_THAT(steep, WithinRel(std::log(1e-3) / 20, 1e-14));
    /// a discontinuous function: the bracket closes on the jump
    double jump = math::findRoot([](double x) { return x < 0.3 ? -1.0 : 1.0; }, 0, 1);
    CHECK_THAT(jump, WithinAbs(0.3, 1e-15));
    CHECK_THROWS_AS(math::findRoot([](double x) { return x * x + 1; }, -1, 1), std::invalid_argument);
}

TEST_CASE("Multivariate minimization", "[numerics]")
{
    auto rosenbrock = [](std::span<const double> x) {
        return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
    };
    auto simplex = math::nelderMead(rosenbrock, {-1.2, 1.0});
    CHECK(simplex.converged);
    CHECK_THAT(simplex.x[0], WithinAbs(1.0, 1e-6));
    CHECK_THAT(simplex.x[1], WithinAbs(1.0, 1e-6));

    auto quasiNewton = math::bfgs(rosenbrock, {-1.2, 1.0});
    CHECK(quasiNewton.converged);
    CHECK_THAT(quasiNewton.x[0], WithinAbs(1.0, 1e-5));
    CHECK_THAT(quasiNewton.x[1], WithinAbs(1.0, 1e-5));

    auto quadratic = [](std::span<const double> x) {
        double sum = 0;
        for (size_t i = 0; i < x.size(); ++i)
            sum += (i + 1) * (x[i] - i) * (x[i] - i);
        return sum;
    };
    auto result = math::bfgs(quadratic, std::vector<double>(5, 10.0));
    for (size_t i = 0; i < 5; ++i)
        CHECK_THAT(result.x[i], WithinAbs(static_cast<double>(i), 1e-6));
}

TEST_CASE("Numerical Hessian recovers quadratic curvature", "[numerics]")
{
    auto quadratic = [](std::span<const double> x) {
        return 2 * x[0] * x[0] + 3 * x[0] * x[1] + 5 * x[1] * x[1] + 7 * x[0] - 11 * x[1];
    };
    auto hessian = math::numericalHessian(quadratic, std::vector<double>{1.25, -0.75});
    REQUIRE(hessian.size() == 2);
    REQUIRE(hessian[0].size() == 2);
    CHECK_THAT(hessian[0][0], WithinAbs(4.0, 1e-6));
    CHECK_THAT(hessian[0][1], WithinAbs(3.0, 1e-6));
    CHECK_THAT(hessian[1][0], WithinAbs(3.0, 1e-6));
    CHECK_THAT(hessian[1][1], WithinAbs(10.0, 1e-6));
}

TEST_CASE("Exponential-family measures support scalar and pair parameters", "[numerics]")
{
    NormalRand<double> normal(1.5, 2.0);
    CHECK_THAT(normal.LogProbabilityMeasure(0.25), WithinRel(normal.logf(0.25), 1e-15));

    constexpr double otherMean = -0.5;
    constexpr double otherVariance = 3.0;
    const double expectedDivergence =
        0.5 * (std::log(otherVariance / normal.Variance()) +
               (normal.Variance() + std::pow(normal.Mean() - otherMean, 2)) / otherVariance - 1);
    CHECK_THAT(normal.KullbackLeiblerDivergence({otherMean, otherVariance}), WithinRel(expectedDivergence, 1e-15));

    PoissonRand<int> poisson(3.5);
    CHECK_THAT(poisson.LogProbabilityMeasure(2), WithinRel(poisson.logP(2), 1e-15));
}
