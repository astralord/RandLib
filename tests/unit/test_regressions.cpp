// Regression tests for specific bugs. IDs B1-B10 are those of docs/MODERNIZATION_PLAN.md,
// the others are listed in CHANGELOG.md under 1.2.0. Expected values come from mpmath
// (tools/refdata) unless they follow from a closed form written in the test.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "RandLib.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
const double PI = 3.14159265358979323846;
}

TEST_CASE("B1: normal samples do not depend on other instances", "[regression]")
{
    /// The ziggurat used to cache a leftover exponential in a static thread_local variable
    NormalRand<double> a(0, 1), b(0, 1), other(0, 1);
    a.Reseed(42);
    b.Reseed(42);
    std::vector<double> first(100000), noise(1000), second(100000);
    a.Sample(first);
    other.Sample(noise);
    b.Sample(second);
    CHECK(first == second);
}

TEST_CASE("B5: the selected generator yields uniforms in (0, 1)", "[regression]")
{
    RandGenerator generator;
    generator.Reseed(5);
    double smallest = 1, largest = 0;
    for (int i = 0; i < 100000; ++i) {
        double u = UniformRand<double>::StandardVariate(generator);
        smallest = std::min(smallest, u);
        largest = std::max(largest, u);
    }
    CHECK(smallest > 0.0);
    CHECK(largest < 1.0);
    CHECK(largest > 0.999);
}

TEST_CASE("B7: stable laws with exponent close to 1 can be constructed and sampled", "[regression]")
{
    StableRand<double> nearOne(0.995, 0.5, 1, 0), smallSkew(1, 0.005, 1, 0);
    nearOne.Reseed(7);
    smallSkew.Reseed(7);
    std::vector<double> sample(1000);
    nearOne.Sample(sample);
    CHECK(std::all_of(sample.begin(), sample.end(), [](double x) { return std::isfinite(x); }));
    smallSkew.Sample(sample);
    CHECK(std::all_of(sample.begin(), sample.end(), [](double x) { return std::isfinite(x); }));
    CHECK(std::abs(nearOne.CF(1.0)) <= 1.0);
    CHECK_THROWS_AS(nearOne.f(0.0), randlib::not_implemented);
    CHECK_THROWS_AS(nearOne.F(0.0), randlib::not_implemented);
    CHECK_THROWS_AS(smallSkew.Quantile(0.3), randlib::not_implemented);
}

TEST_CASE("B8: totally skewed stable distribution function has no floor", "[regression]")
{
    StableRand<double> X(0.6, 1, 1, 0);
    CHECK(X.F(0.0) == 0.0);
    CHECK(X.F(0.01) < 1e-30);
    CHECK_THAT(X.F(0.2542), WithinRel(8.5586231489170863891e-4, 1e-9));
}

TEST_CASE("B9: normal samples are symmetric", "[regression]")
{
    NormalRand<double> X(0, 1);
    X.Reseed(9);
    std::vector<double> sample(1000000);
    X.Sample(sample);
    double positive = static_cast<double>(std::count_if(sample.begin(), sample.end(), [](double x) { return x > 0; }));
    /// 5 standard deviations of a Binomial(10^6, 1/2) count
    CHECK(std::fabs(positive - 500000.0) < 2500.0);
}

TEST_CASE("B11: quantiles are deterministic", "[regression]")
{
    /// The generic quantile used a sampled initial guess from a shared static vector
    BetaPrimeRand<double> a(2, 3), b(2, 3);
    a.Reseed(1);
    b.Reseed(2);
    CHECK(a.Quantile(0.3) == b.Quantile(0.3));
    CHECK(a.Quantile(0.3) == a.Quantile(0.3));
}

TEST_CASE("B12: Bernoulli bits come from the generator that is passed", "[regression]")
{
    RandGenerator first, second;
    first.Reseed(12);
    second.Reseed(12);
    std::vector<int> a, b;
    for (int i = 0; i < 1000; ++i) {
        a.push_back(BernoulliRand::StandardVariate(first));
        b.push_back(BernoulliRand::StandardVariate(second));
    }
    CHECK(a == b);
}

TEST_CASE("B13: logarithm of pi has full precision", "[regression]")
{
    CHECK_THAT(static_cast<double>(randlib::consts::lnpi), WithinRel(std::log(PI), 1e-16));
    CauchyRand<double> X(0, 1);
    CHECK_THAT(X.logf(0.0), WithinRel(-std::log(PI), 1e-15));
}

TEST_CASE("B14: normal raw moments", "[regression]")
{
    NormalRand<double> X(2, 0.09);
    double mu = 2, var = 0.09;
    CHECK_THAT(static_cast<double>(X.ThirdMoment()), WithinRel(mu * mu * mu + 3 * mu * var, 1e-14));
    CHECK_THAT(static_cast<double>(X.FourthMoment()), WithinRel(std::pow(mu, 4) + 6 * mu * mu * var + 3 * var * var, 1e-14));
}

TEST_CASE("B15: classes that could not be used", "[regression]")
{
    /// BaldingNicholsRand hid F(x) behind a data member, CenteredAsymmetricLaplaceRand was abstract,
    /// WrappedExponentialRand was never instantiated and base destructors were protected
    BaldingNicholsRand<double> balding(0.3, 0.6);
    CHECK(balding.F(0.5) > 0.0);
    CenteredAsymmetricLaplaceRand<double> laplace(1, 2);
    CHECK(laplace.f(0.0) > 0.0);
    WrappedExponentialRand<double> wrapped(1);
    CHECK(wrapped.F(1.0) > 0.0);
    std::unique_ptr<UnivariateDistribution<double>> owned = std::make_unique<NormalRand<double>>(0, 1);
    CHECK(owned->Mean() == 0.0);
}

TEST_CASE("B16: wrapped exponential circular moments and characteristic function", "[regression]")
{
    double lambda = 0.7;
    WrappedExponentialRand<double> X(lambda);
    CHECK_THAT(static_cast<double>(X.CircularVariance()), WithinRel(1 - lambda / std::sqrt(1 + lambda * lambda), 1e-14));
    std::complex<double> expected = lambda / std::complex<double>(lambda, -1.0);
    CHECK_THAT(std::abs(X.CF(1.0) - expected), WithinAbs(0, 1e-14));
    CHECK_THROWS_AS(X.SetRate(-1), std::invalid_argument);
}

TEST_CASE("B17: incomplete gamma ratios keep full precision for large shapes", "[regression]")
{
    /// Poisson(40) and Gamma(50, 3) used to be off by up to 1e-3
    PoissonRand<int> P(40);
    CHECK_THAT(P.F(40), WithinRel(0.54191817836253703985, 1e-13));
    GammaRand<double> G(50, 3);
    CHECK_THAT(G.F(15.0), WithinRel(0.24680203440017027271, 1e-12));
}

TEST_CASE("B18: logarithm of the beta function does not underflow", "[regression]")
{
    CHECK_THAT(static_cast<double>(RandMath::logBeta(600, 600)), WithinRel(-833.70936104274857839, 1e-15));
    BetaRand<double> X(600, 600);
    CHECK(std::isfinite(X.f(0.5)));
}

TEST_CASE("B19: right-triangular distributions", "[regression]")
{
    TriangularRand<double> X(0, 0, 1);
    CHECK(X.f(0.0) == 2.0);
    CHECK_THAT(X.F(0.5), WithinRel(0.75, 1e-15));
    CHECK(std::isfinite(X.CF(1.0).real()));
}

TEST_CASE("B20: modified Bessel function of negative and half-integer order", "[regression]")
{
    /// I_{1/2}(x) = sqrt(2 / (pi x)) sinh x and I_{-1/2}(x) = sqrt(2 / (pi x)) cosh x
    double x = 2.0;
    CHECK_THAT(static_cast<double>(RandMath::logBesselI(0.5, x)), WithinRel(std::log(std::sinh(x) * std::sqrt(2 / (PI * x))), 1e-14));
    CHECK_THAT(static_cast<double>(RandMath::logBesselI(-0.5, x)), WithinRel(std::log(std::cosh(x) * std::sqrt(2 / (PI * x))), 1e-14));
}

TEST_CASE("B21: stable characteristic function and location for exponent 1", "[regression]")
{
    /// S1 parameterization: phi(t) = exp(-(gamma t)^a (1 - i b tan(pi a / 2)) + i mu t) for t > 0
    StableRand<double> X(1.5, 0.5, 1, 0);
    std::complex<double> expected = std::exp(std::complex<double>(-1.0, 0.5 * std::tan(PI * 0.75)));
    CHECK_THAT(std::abs(X.CF(1.0) - expected), WithinAbs(0, 1e-14));
    /// alpha = 1 with gamma != 1: the S1 location shift is +(2/pi) b gamma log gamma
    StableRand<double> Y(1, 0.5, 2, 1);
    CHECK_THAT(Y.F(1.0), WithinRel(0.37109122054964755289, 1e-9));
}

TEST_CASE("B22: stable tails", "[regression]")
{
    StableRand<double> unity(1, 0.5, 1, 0);
    CHECK_THAT(unity.f(300.0), WithinRel(5.3592152809976668665e-6, 1e-9));
    StableRand<double> skewed(1.3, 1, 1, 0);
    CHECK_THAT(skewed.f(134.3), WithinRel(8.3768206528387327408e-6, 1e-9));
    CHECK_THAT(skewed.S(638542.0), WithinRel(1.445556098e-8, 1e-8));
}

TEST_CASE("B23: survival functions in the far tails", "[regression]")
{
    CauchyRand<double> cauchy(0, 1);
    CHECK_THAT(cauchy.F(-3.183098861833411e29), WithinRel(1e-30, 1e-9));
    BetaPrimeRand<double> betaPrime(2, 3);
    CHECK(betaPrime.S(3.4199518933533938e33) > 0.0);
    StudentTRand<double> student(1.5, 1, 2);
    CHECK(student.F(-1.0438938854689273e200) > 0.0);
    CHECK(std::isfinite(student.logf(-1.0438938854689273e200)));
    LogisticRand<double> logistic(0, 1);
    CHECK(logistic.f(-690.77552789821368) > 0.0);
}

TEST_CASE("B24: discrete distribution functions and medians", "[regression]")
{
    LogarithmicRand<int> logarithmic(0.9);
    CHECK_THAT(logarithmic.F(3) - logarithmic.F(2), WithinRel(logarithmic.P(3), 1e-13));
    BinomialRand<int> binomial(5, 0.9);
    CHECK(binomial.Median() == 5);
    NegativeBinomialRand<int> negative(2.5, 0.4);
    CHECK(negative.S(-1) == 1.0);
    BernoulliRand bernoulli(0.3);
    CHECK(bernoulli.logP(2) == -INFINITY);
    CategoricalRand<int> categorical({0.2, 0.3, 0.5});
    CHECK(categorical.Quantile1m(0.5) == 1);
}

TEST_CASE("B25: invalid parameters are not silently replaced", "[regression]")
{
    CHECK_THROWS_AS(LogNormalRand<double>(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(BaldingNicholsRand<double>(1.5, 0.5), std::invalid_argument);
    CHECK_THROWS_AS(NegativeBinomialRand<int>(-1, 0.5), std::invalid_argument);
}

TEST_CASE("B26: composite distributions reseed every generator they use", "[regression]")
{
    PlanckRand<double> a(3, 2), b(3, 2);
    a.Reseed(26);
    b.Reseed(26);
    CHECK(a.Variate() == b.Variate());
    LaplaceRand<double> c(0, 1), d(0, 1);
    c.Reseed(26);
    d.Reseed(26);
    CHECK(c.Variate() == d.Variate());
    CHECK(RandEngine::DeriveSeed(42, 1) != 43UL);
}

TEST_CASE("B27: heavy-tailed discrete characteristic functions terminate", "[regression]")
{
    ZetaRand<int> zeta(2);
    auto start = std::chrono::steady_clock::now();
    std::complex<double> phi = zeta.CF(0.5);
    CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(10));
    /// Li_2(e^{0.5i}) / zeta(2)
    CHECK_THAT(phi.real(), WithinAbs(0.56053061459019065698, 1e-12));
    CHECK_THAT(phi.imag(), WithinAbs(0.51571178128074520922, 1e-12));
    CHECK(zeta.Quantile1m(1e-10) > 0);

    /// the tail is summed with Euler's transformation instead of term by term: Li_s(e^{it}) / zeta(s)
    ZetaRand<int> heavy(1.1), light(3.5);
    start = std::chrono::steady_clock::now();
    std::complex<double> heavyPhi = heavy.CF(1.0), lightPhi = light.CF(0.1);
    CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(1));
    CHECK_THAT(heavyPhi.real(), WithinAbs(0.0076240340145568357863, 1e-14));
    CHECK_THAT(heavyPhi.imag(), WithinAbs(0.10101969240488333018, 1e-14));
    CHECK_THAT(lightPhi.real(), WithinAbs(0.99028256032171156432, 1e-14));
    CHECK_THAT(lightPhi.imag(), WithinAbs(0.1173998095713813435, 1e-14));
}

TEST_CASE("B28: Cantor characteristic function depends on t", "[regression]")
{
    CantorRand X;
    double t = 3.0, product = 1.0;
    for (int k = 1; k < 40; ++k)
        product *= std::cos(t / std::pow(3.0, k));
    CHECK_THAT(X.CF(t).real(), WithinAbs(std::cos(t / 2) * product, 1e-12));
}

TEST_CASE("B29: geometric stable distribution function for exponent 1", "[regression]")
{
    /// The mixture integral missed the alpha = 1 shift (2/pi) b gamma w log(gamma w)
    GeometricStableRand<double> X(1, 0.5, 1, 0);
    X.Reseed(29);
    std::vector<double> sample(20000);
    X.Sample(sample);
    double x = 1.0;
    double empirical = static_cast<double>(std::count_if(sample.begin(), sample.end(), [x](double v) { return v <= x; })) / 20000.0;
    CHECK(std::fabs(empirical - X.F(x)) < 0.02);
}

TEST_CASE("B30: Marchenko-Pastur distribution function next to the lower edge", "[regression]")
{
    /// For x below ~1e-16 the square root in the formula vanished and F(x) came out as 1/4
    MarchenkoPasturRand<double> X(1, 1);
    CHECK(X.F(1e-17) < 1e-8);
    CHECK_THAT(X.F(1e-10), WithinRel(2 / PI * 1e-5, 1e-4));
}

TEST_CASE("B31: Rademacher mode is deterministic", "[regression]")
{
    RademacherRand X;
    for (int i = 0; i < 20; ++i)
        CHECK(X.Mode() == -1);
}

TEST_CASE("B32: Poisson and binomial probabilities do not depend on the accuracy of lgamma", "[regression]")
{
    /// k log(λ) - λ - log(k!) cancelled about four digits; the saddle-point form keeps full precision
    CHECK_THAT(PoissonRand<int>(1000).P(977), WithinRel(0.009776061021570442222212799, 1e-14));
    CHECK_THAT(BinomialRand<int>(1000, 0.5).P(480), WithinRel(0.01134098237118772746599741, 1e-14));
    CHECK_THAT(BinomialRand<int>(100, 0.02).P(7), WithinRel(0.003130114101451260277632229, 1e-14));
}

TEST_CASE("B33: raised cosine characteristic function at st = pi", "[regression]")
{
    /// (1 + cos x) / (2 pi) on [-pi, pi]: CF(1) = 1/2, the limit of a 0 / 0 expression
    CHECK_THAT(RaisedCosineRand<double>(0, PI).CF(1).real(), WithinAbs(0.5, 1e-15));
    CHECK_THAT(RaabGreenRand<double>().CF(1).real(), WithinAbs(0.5, 1e-15));
    CHECK_THAT(RaisedCosineRand<double>(0, PI).CF(1 + 1e-9).real(), WithinAbs(0.5, 1e-8));
}

TEST_CASE("B34: discrete quantiles at an exact probability level", "[regression]")
{
    /// F(2) = 1/2 exactly for NB(3, 1/2); a last-bit error in F must not move the median
    CHECK(PascalRand<int>(3, 0.5).Quantile(0.5) == 2);
    CHECK(PascalRand<int>(3, 0.5).Median() == 2);
}

TEST_CASE("B35: zeta probabilities in the far tail terminate", "[regression]")
{
    ZetaRand<int> X(3.5);
    /// mpmath Hurwitz-zeta references at 80 decimal digits
    CHECK_THAT(X.S(100001), WithinRel(1.122593069926359984040757874880502e-13, 1e-15));
    CHECK_THAT(X.S(1000000000), WithinRel(1.122635166277711578412826611568585e-23, 1e-15));
}

TEST_CASE("B36: uniform likelihoods use a negative sample-count exponent", "[regression]")
{
    /// The exact likelihood is p(x)^n; negating the unsigned sample size wrapped to a huge positive exponent
    const std::vector<double> continuousSample = {0.25, 0.5, 1.5};
    CHECK(UniformRand<double>(0, 2).LikelihoodFunction(continuousSample) == 0.125);

    const std::vector<int> discreteSample = {0, 1, 0};
    CHECK(UniformDiscreteRand<int>(0, 1).LikelihoodFunction(discreteSample) == 0.125);
}

TEST_CASE("B37: double factorial uses every second factor", "[regression]")
{
    CHECK(RandMath::doubleFactorial(-1) == 0.0);
    CHECK(RandMath::doubleFactorial(0) == 1.0);
    CHECK(RandMath::doubleFactorial(1) == 1.0);
    CHECK(RandMath::doubleFactorial(8) == 384.0);
    CHECK(RandMath::doubleFactorial(9) == 945.0);
}

TEST_CASE("B38: geometric waiting times do not overflow signed integers", "[regression]")
{
    const double probability = std::numeric_limits<double>::min();
    RandGenerator generator;
    generator.Reseed(20231010);

    CHECK(GeometricRand<int>::Variate(probability, generator) == std::numeric_limits<int>::max());

    GeometricRand<int> geometric(probability);
    geometric.Reseed(20231010);
    CHECK(geometric.Variate() == std::numeric_limits<int>::max());

    PascalRand<int> negativeBinomial(2, probability);
    negativeBinomial.Reseed(20231010);
    CHECK(negativeBinomial.Variate() == std::numeric_limits<int>::max());

    generator.Reseed(20231010);
    CHECK(BinomialRand<int>::Variate(20, probability, generator) == 0);
}

TEST_CASE("B39: negative-binomial fitting stays in the positive shape domain", "[regression]")
{
    double invalidRoot = 1.0;
    int evaluations = 0;
    CHECK_FALSE(RandMath::findRootNewtonFirstOrder<double>(
        [&evaluations] (double) {
            ++evaluations;
            return DoublePair(1.0, std::numeric_limits<double>::quiet_NaN());
        },
        invalidRoot));
    CHECK(evaluations == 1);

    NegativeBinomialRand<int> fitted;
    fitted.Fit({0, 1, 3});
    CHECK_THAT(fitted.GetNumber(), WithinRel(5.2122598849254577174, 1e-10));
    CHECK_THAT(fitted.GetProbability(), WithinRel(0.79630061189656137372, 1e-10));
}

TEST_CASE("B40: Rayleigh unbiased scale fitting stays in the log domain", "[regression]")
{
    /// If R = sqrt(sum(X_i^2) / (2n)), then E[R] = sigma Gamma(n + 1/2) / (sqrt(n) Gamma(n)).
    for (size_t n : {size_t{15}, size_t{16}, size_t{101}}) {
        const std::vector<double> sample(n, 1.0);
        const double mleScaleSquared = 0.5;
        const double expected =
            std::exp(0.5 * std::log(n * mleScaleSquared) + std::lgamma(n) - std::lgamma(n + 0.5));

        RayleighRand<double> fitted;
        fitted.Fit(sample, true);

        INFO("sample size " << n);
        CHECK_THAT(fitted.GetScale(), WithinRel(expected, 1e-12));
    }
}

TEST_CASE("B41: asymmetric-Laplace entropy includes the exponential contribution", "[regression]")
{
    /// Integrating -log(f) gives 1 + log(gamma * (kappa + 1 / kappa)).
    CHECK_THAT(static_cast<double>(LaplaceRand<double>(0, 1).Entropy()), WithinRel(1 + std::log(2.0), 1e-15));
    CHECK_THAT(static_cast<double>(LaplaceRand<double>(4, 2).Entropy()), WithinRel(1 + std::log(4.0), 1e-15));
    CHECK_THAT(static_cast<double>(AsymmetricLaplaceRand<double>(7, 3, 2).Entropy()),
               WithinRel(1 + std::log(7.5), 1e-15));
}

TEST_CASE("B42: zero noncentrality is the central chi-squared limit", "[regression]")
{
    const NoncentralChiSquaredRand<double> central(4, 0);
    CHECK_THAT(central.f(2), WithinRel(0.5 / std::exp(1.0), 1e-14));
    CHECK_THAT(central.F(2), WithinRel(1 - 2 / std::exp(1.0), 1e-14));
    CHECK_THAT(central.S(2), WithinRel(2 / std::exp(1.0), 1e-14));
    CHECK_THAT(central.CF(0.5).real(), WithinAbs(0.0, 1e-15));
    CHECK_THAT(central.CF(0.5).imag(), WithinAbs(0.5, 1e-15));
    CHECK(central.Mean() == 4);
    CHECK(central.Variance() == 8);
    CHECK(central.Mode() == 2);

    RandGenerator actualGenerator, expectedGenerator;
    actualGenerator.Reseed(42);
    expectedGenerator.Reseed(42);
    const double actual = NoncentralChiSquaredRand<double>::Variate(0.5, 0, actualGenerator);
    const double expected = 2 * GammaDistribution<double>::StandardVariate(0.25, expectedGenerator);
    CHECK(actual == expected);
}

TEST_CASE("B43: Bernoulli endpoint entropy uses the zero-log-zero limit", "[regression]")
{
    CHECK(BernoulliRand(0).Entropy() == 0);
    CHECK(BernoulliRand(1).Entropy() == 0);
    CHECK_THAT(BernoulliRand(0.5).Entropy(), WithinRel(std::log(2.0), 1e-15));
}

TEST_CASE("B44: binomial probability boundaries preserve construction and sampling", "[regression]")
{
    const BinomialRand<int> certain(7, 1);
    CHECK(certain.Mode() == 7);
    CHECK(certain.Variate() == 7);
    std::vector<int> sample(3);
    certain.Sample(sample);
    CHECK(sample == std::vector<int>{7, 7, 7});

    const BinomialRand<int> impossible(7, 0);
    CHECK(impossible.Mode() == 0);
    CHECK(impossible.Variate() == 0);
    CHECK(BernoulliRand(1).Mode() == 1);
    CHECK(BinomialRand<int>(5, 0.5).Mode() == 3);

    BinomialRand<int> bulk(10000000, 1 - 9e-7), scalar(10000000, 1 - 9e-7);
    bulk.Reseed(42);
    scalar.Reseed(42);
    std::vector<int> actual(4), expected(4);
    bulk.Sample(actual);
    for (int &value : expected)
        value = scalar.Variate();
    CHECK(actual == expected);
}
