// Expected values come from mpmath (tools/refdata) unless they follow from a closed form written
// in the test.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

using namespace randlib;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
const double PI = 3.14159265358979323846;
}

TEST_CASE("normal samples do not depend on other instances", "[regression]")
{
    /// Separate distribution objects must not share cached sampler state.
    NormalRand<double> a(0, 1), b(0, 1), other(0, 1);
    a.Reseed(42);
    b.Reseed(42);
    std::vector<double> first(100000), noise(1000), second(100000);
    a.Sample(first);
    other.Sample(noise);
    b.Sample(second);
    CHECK(first == second);
}

TEST_CASE("64-bit engines yield uniforms in (0, 1)", "[regression]")
{
    DefaultEngine generator(5);
    double smallest = 1, largest = 0;
    for (int i = 0; i < 100000; ++i) {
        double u = detail::uniform01<double>(generator);
        smallest = std::min(smallest, u);
        largest = std::max(largest, u);
    }
    CHECK(smallest > 0.0);
    CHECK(largest < 1.0);
    CHECK(largest > 0.999);
}

TEST_CASE("stable laws with exponent close to 1 can be constructed and sampled", "[regression]")
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
    CHECK(std::isfinite(nearOne.f(0.0)));
    CHECK(nearOne.F(0.0) > 0.0);
    CHECK(std::isfinite(smallSkew.Quantile(0.3)));
}

TEST_CASE("totally skewed stable distribution function has no floor", "[regression]")
{
    StableRand<double> X(0.6, 1, 1, 0);
    CHECK(X.F(0.0) == 0.0);
    CHECK(X.F(0.01) < 1e-30);
    CHECK_THAT(X.F(0.2542), WithinRel(8.5586231489170863891e-4, 1e-9));
}

TEST_CASE("normal samples are symmetric", "[regression]")
{
    NormalRand<double> X(0, 1);
    X.Reseed(9);
    std::vector<double> sample(1000000);
    X.Sample(sample);
    double positive = static_cast<double>(std::count_if(sample.begin(), sample.end(), [](double x) { return x > 0; }));
    /// 5 standard deviations of a Binomial(10^6, 1/2) count
    CHECK(std::fabs(positive - 500000.0) < 2500.0);
}

TEST_CASE("quantiles are deterministic", "[regression]")
{
    /// Quantile evaluation must be deterministic and free of shared sampler state.
    BetaPrimeRand<double> a(2, 3), b(2, 3);
    a.Reseed(1);
    b.Reseed(2);
    CHECK(a.Quantile(0.3) == b.Quantile(0.3));
    CHECK(a.Quantile(0.3) == a.Quantile(0.3));
}

TEST_CASE("Bernoulli bits come from the generator that is passed", "[regression]")
{
    Pcg32 first(12), second(12);
    std::vector<int> a, b;
    for (int i = 0; i < 1000; ++i) {
        a.push_back(BernoulliRand::StandardVariate(first));
        b.push_back(BernoulliRand::StandardVariate(second));
    }
    CHECK(a == b);
}

TEST_CASE("logarithm of pi has full precision", "[regression]")
{
    CHECK_THAT(static_cast<double>(randlib::consts::lnpi), WithinRel(std::log(PI), 1e-16));
    CauchyRand<double> X(0, 1);
    CHECK_THAT(X.logf(0.0), WithinRel(-std::log(PI), 1e-15));
}

TEST_CASE("normal raw moments", "[regression]")
{
    NormalRand<double> X(2, 0.09);
    double mu = 2, var = 0.09;
    CHECK_THAT(static_cast<double>(X.ThirdMoment()), WithinRel(mu * mu * mu + 3 * mu * var, 1e-14));
    CHECK_THAT(static_cast<double>(X.FourthMoment()),
               WithinRel(std::pow(mu, 4) + 6 * mu * mu * var + 3 * var * var, 1e-14));
}

TEST_CASE("classes that could not be used", "[regression]")
{
    /// Public concrete distributions support probability evaluation and polymorphic destruction.
    BaldingNicholsRand<double> balding(0.3, 0.6);
    CHECK(balding.F(0.5) > 0.0);
    CenteredAsymmetricLaplaceRand<double> laplace(1, 2);
    CHECK(laplace.f(0.0) > 0.0);
    WrappedExponentialRand<double> wrapped(1);
    CHECK(wrapped.F(1.0) > 0.0);
    std::unique_ptr<UnivariateDistribution<double>> owned = std::make_unique<NormalRand<double>>(0, 1);
    CHECK(owned->Mean() == 0.0);
}

TEST_CASE("wrapped exponential circular moments and characteristic function", "[regression]")
{
    double lambda = 0.7;
    WrappedExponentialRand<double> X(lambda);
    CHECK_THAT(static_cast<double>(X.CircularVariance()),
               WithinRel(1 - lambda / std::sqrt(1 + lambda * lambda), 1e-14));
    std::complex<double> expected = lambda / std::complex<double>(lambda, -1.0);
    CHECK_THAT(std::abs(X.CF(1.0) - expected), WithinAbs(0, 1e-14));
    CHECK_THROWS_AS(X.SetRate(-1), std::invalid_argument);
}

TEST_CASE("incomplete gamma ratios keep full precision for large shapes", "[regression]")
{
    /// Large-shape incomplete-gamma ratios retain full precision.
    PoissonRand<int> P(40);
    CHECK_THAT(P.F(40), WithinRel(0.54191817836253703985, 1e-13));
    GammaRand<double> G(50, 3);
    CHECK_THAT(G.F(15.0), WithinRel(0.24680203440017027271, 1e-12));
}

TEST_CASE("logarithm of the beta function does not underflow", "[regression]")
{
    CHECK_THAT(static_cast<double>(math::logBeta(600, 600)), WithinRel(-833.70936104274857839, 1e-15));
    BetaRand<double> X(600, 600);
    CHECK(std::isfinite(X.f(0.5)));
}

TEST_CASE("right-triangular distributions", "[regression]")
{
    TriangularRand<double> X(0, 0, 1);
    CHECK(X.f(0.0) == 2.0);
    CHECK_THAT(X.F(0.5), WithinRel(0.75, 1e-15));
    CHECK(std::isfinite(X.CF(1.0).real()));
}

TEST_CASE("modified Bessel function of negative and half-integer order", "[regression]")
{
    /// I_{1/2}(x) = sqrt(2 / (pi x)) sinh x and I_{-1/2}(x) = sqrt(2 / (pi x)) cosh x
    double x = 2.0;
    CHECK_THAT(static_cast<double>(math::logBesselI(0.5, x)),
               WithinRel(std::log(std::sinh(x) * std::sqrt(2 / (PI * x))), 1e-14));
    CHECK_THAT(static_cast<double>(math::logBesselI(-0.5, x)),
               WithinRel(std::log(std::cosh(x) * std::sqrt(2 / (PI * x))), 1e-14));
}

TEST_CASE("stable characteristic function and location for exponent 1", "[regression]")
{
    /// S1 parameterization: phi(t) = exp(-(gamma t)^a (1 - i b tan(pi a / 2)) + i mu t) for t > 0
    StableRand<double> X(1.5, 0.5, 1, 0);
    std::complex<double> expected = std::exp(std::complex<double>(-1.0, 0.5 * std::tan(PI * 0.75)));
    CHECK_THAT(std::abs(X.CF(1.0) - expected), WithinAbs(0, 1e-14));
    /// alpha = 1 with gamma != 1: the S1 location shift is +(2/pi) b gamma log gamma
    StableRand<double> Y(1, 0.5, 2, 1);
    CHECK_THAT(Y.F(1.0), WithinRel(0.37109122054964755289, 1e-9));
}

TEST_CASE("stable tails", "[regression]")
{
    StableRand<double> unity(1, 0.5, 1, 0);
    CHECK_THAT(unity.f(300.0), WithinRel(5.3592152809976668665e-6, 1e-9));
    StableRand<double> skewed(1.3, 1, 1, 0);
    CHECK_THAT(skewed.f(134.3), WithinRel(8.3768206528387327408e-6, 1e-9));
    CHECK_THAT(skewed.S(638542.0), WithinRel(1.445556098e-8, 1e-8));
}

TEST_CASE("survival functions in the far tails", "[regression]")
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

TEST_CASE("discrete distribution functions and medians", "[regression]")
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

TEST_CASE("categorical guide-table inversion equals sequential search", "[regression]")
{
    std::vector<double> probabilities(1000);
    for (std::size_t i = 0; i < probabilities.size(); ++i)
        probabilities[i] = 1.0 + static_cast<double>((i * 7919) % 101);
    probabilities.front() = probabilities.back() = 0;
    const double total = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    for (double &probability : probabilities)
        probability /= total;
    CategoricalRand<int> categorical{std::vector<double>(probabilities)};
    auto sequential = [&probabilities](double p) {
        double cumulative = 0;
        int k = 0;
        while ((cumulative += probabilities[k]) < p && k + 1 < static_cast<int>(probabilities.size()))
            ++k;
        return k;
    };
    std::mt19937_64 engine(20231010);
    std::uniform_real_distribution<double> uniform;
    for (int i = 0; i < 20000; ++i) {
        const double p = uniform(engine);
        CHECK(categorical.Quantile(p) == sequential(p));
    }
    for (double p : {1e-300, 0.5, 1 - 1e-16})
        CHECK(categorical.Quantile(p) == sequential(p));
}

TEST_CASE("invalid parameters are not silently replaced", "[regression]")
{
    CHECK_THROWS_AS(LogNormalRand<double>(0, -1), std::invalid_argument);
    CHECK_THROWS_AS(BaldingNicholsRand<double>(1.5, 0.5), std::invalid_argument);
    CHECK_THROWS_AS(NegativeBinomialRand<int>(-1, 0.5), std::invalid_argument);
}

TEST_CASE("composite distributions reseed every generator they use", "[regression]")
{
    PlanckRand<double> a(3, 2), b(3, 2);
    a.Reseed(26);
    b.Reseed(26);
    CHECK(a.Variate() == b.Variate());
    LaplaceRand<double> c(0, 1), d(0, 1);
    c.Reseed(26);
    d.Reseed(26);
    CHECK(c.Variate() == d.Variate());
}

TEST_CASE("heavy-tailed discrete characteristic functions terminate", "[regression]")
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

TEST_CASE("Cantor characteristic function depends on t", "[regression]")
{
    CantorRand X;
    double t = 3.0, product = 1.0;
    for (int k = 1; k < 40; ++k)
        product *= std::cos(t / std::pow(3.0, k));
    CHECK_THAT(X.CF(t).real(), WithinAbs(std::cos(t / 2) * product, 1e-12));
}

TEST_CASE("geometric stable distribution function for exponent 1", "[regression]")
{
    /// The alpha = 1 mixture includes the shift (2/pi) b gamma w log(gamma w).
    GeometricStableRand<double> X(1, 0.5, 1, 0);
    X.Reseed(29);
    std::vector<double> sample(20000);
    X.Sample(sample);
    double x = 1.0;
    double empirical =
        static_cast<double>(std::count_if(sample.begin(), sample.end(), [x](double v) { return v <= x; })) / 20000.0;
    CHECK(std::fabs(empirical - X.F(x)) < 0.02);

    // Direct characteristic-function asymptotics evaluated with mpmath at 40 decimal digits.
    CHECK_THAT(X.f(-119), WithinAbs(0.00001079183336235705390278076, 1e-9));
    CHECK_THAT(X.F(-119), WithinAbs(0.0013070736358547730840764597, 1e-9));
    CHECK_THAT(X.f(570.9), WithinAbs(0.000001482821763075677554750345, 1e-9));
    CHECK_THAT(X.F(570.9), WithinAbs(0.9991581026900738992757455, 1e-9));
}

TEST_CASE("Marchenko-Pastur distribution function next to the lower edge", "[regression]")
{
    /// The lower-edge expansion remains accurate after direct square-root terms lose resolution.
    MarchenkoPasturRand<double> X(1, 1);
    CHECK(X.F(1e-17) < 1e-8);
    CHECK_THAT(X.F(1e-10), WithinRel(2 / PI * 1e-5, 1e-4));
}

TEST_CASE("Rademacher mode is deterministic", "[regression]")
{
    RademacherRand X;
    for (int i = 0; i < 20; ++i)
        CHECK(X.Mode() == -1);
}

TEST_CASE("Poisson and binomial probabilities do not depend on the accuracy of lgamma", "[regression]")
{
    /// k log(λ) - λ - log(k!) cancelled about four digits; the saddle-point form keeps full precision
    CHECK_THAT(PoissonRand<int>(1000).P(977), WithinRel(0.009776061021570442222212799, 1e-14));
    CHECK_THAT(BinomialRand<int>(1000, 0.5).P(480), WithinRel(0.01134098237118772746599741, 1e-14));
    CHECK_THAT(BinomialRand<int>(100, 0.02).P(7), WithinRel(0.003130114101451260277632229, 1e-14));
}

TEST_CASE("Cantor quantiles are exact", "[regression]")
{
    /// the inverse of the Cantor function maps the binary digits of p to ternary digits 0 and 2
    auto exact = [](double p) {
        double x = 0, power = 1;
        for (int k = 0; k < 60; ++k) {
            p *= 2;
            power /= 3;
            if (p >= 1) {
                x += 2 * power;
                p -= 1;
            }
        }
        return x;
    };
    CantorRand X;
    for (double p : {1e-6, 0.001, 0.3, 0.7})
        CHECK_THAT(X.Quantile(p), WithinRel(exact(p), 1e-12));
}

TEST_CASE("raised cosine characteristic function at st = pi", "[regression]")
{
    /// (1 + cos x) / (2 pi) on [-pi, pi]: CF(1) = 1/2, the limit of a 0 / 0 expression
    CHECK_THAT(RaisedCosineRand<double>(0, PI).CF(1).real(), WithinAbs(0.5, 1e-15));
    CHECK_THAT(RaabGreenRand<double>().CF(1).real(), WithinAbs(0.5, 1e-15));
    CHECK_THAT(RaisedCosineRand<double>(0, PI).CF(1 + 1e-9).real(), WithinAbs(0.5, 1e-8));
}

TEST_CASE("discrete quantiles at an exact probability level", "[regression]")
{
    /// F(2) = 1/2 exactly for NB(3, 1/2); a last-bit error in F must not move the median
    CHECK(PascalRand<int>(3, 0.5).Quantile(0.5) == 2);
    CHECK(PascalRand<int>(3, 0.5).Median() == 2);
}

TEST_CASE("zeta probabilities in the far tail terminate", "[regression]")
{
    ZetaRand<int> X(3.5);
    /// mpmath Hurwitz-zeta references at 80 decimal digits
    CHECK_THAT(X.S(100001), WithinRel(1.122593069926359984040757874880502e-13, 1e-15));
    CHECK_THAT(X.S(1000000000), WithinRel(1.122635166277711578412826611568585e-23, 1e-15));
}

TEST_CASE("uniform likelihoods use a negative sample-count exponent", "[regression]")
{
    /// The exact likelihood is p(x)^n with a signed negative exponent after taking reciprocals.
    const std::vector<double> continuousSample = {0.25, 0.5, 1.5};
    CHECK(UniformRand<double>(0, 2).LikelihoodFunction(continuousSample) == 0.125);

    const std::vector<int> discreteSample = {0, 1, 0};
    CHECK(UniformDiscreteRand<int>(0, 1).LikelihoodFunction(discreteSample) == 0.125);
}

TEST_CASE("double factorial uses every second factor", "[regression]")
{
    CHECK(math::doubleFactorial(-1) == 0.0);
    CHECK(math::doubleFactorial(0) == 1.0);
    CHECK(math::doubleFactorial(1) == 1.0);
    CHECK(math::doubleFactorial(8) == 384.0);
    CHECK(math::doubleFactorial(9) == 945.0);
}

TEST_CASE("geometric waiting times do not overflow signed integers", "[regression]")
{
    const double probability = std::numeric_limits<double>::min();
    std::mt19937_64 generator(20231010);

    CHECK(GeometricRand<int>::Variate(probability, generator) == std::numeric_limits<int>::max());

    GeometricRand<int> geometric(probability);
    generator.seed(20231010);
    CHECK(geometric(generator) == std::numeric_limits<int>::max());

    generator.seed(20231010);
    PascalRand<int> negativeBinomial(2, probability);
    CHECK(negativeBinomial(generator) == std::numeric_limits<int>::max());

    generator.seed(20231010);
    CHECK(BinomialRand<int>::Variate(20, probability, generator) == 0);
}

TEST_CASE("negative-binomial fitting stays in the positive shape domain", "[regression]")
{
    double invalidRoot = 1.0;
    int evaluations = 0;
    CHECK_FALSE(math::findRootNewtonFirstOrder<double>(
        [&evaluations](double) {
            ++evaluations;
            return DoublePair(1.0, std::numeric_limits<double>::quiet_NaN());
        },
        invalidRoot));
    CHECK(evaluations == 1);

    const std::vector<int> sample = {0, 1, 3};
    auto fitted = fit<PolyaRand<int>>(sample);
    REQUIRE(fitted.converged);
    CHECK_THAT(fitted.distribution.GetNumber(), WithinRel(5.2122598849254577174, 1e-10));
    CHECK_THAT(fitted.distribution.GetProbability(), WithinRel(0.79630061189656137372, 1e-10));
}

TEST_CASE("Rayleigh unbiased scale fitting stays in the log domain", "[regression]")
{
    /// If R = sqrt(sum(X_i^2) / (2n)), then E[R] = sigma Gamma(n + 1/2) / (sqrt(n) Gamma(n)).
    for (size_t n : {size_t{15}, size_t{16}, size_t{101}}) {
        const std::vector<double> sample(n, 1.0);
        const double mleScaleSquared = 0.5;
        const double expected = std::exp(0.5 * std::log(n * mleScaleSquared) + std::lgamma(n) - std::lgamma(n + 0.5));

        const auto fitted = fit<RayleighRand<double>>(sample, Method::UMVU);

        INFO("sample size " << n);
        CHECK_THAT(fitted.distribution.GetScale(), WithinRel(expected, 1e-12));
    }
}

TEST_CASE("asymmetric-Laplace entropy includes the exponential contribution", "[regression]")
{
    /// Integrating -log(f) gives 1 + log(gamma * (kappa + 1 / kappa)).
    CHECK_THAT(static_cast<double>(LaplaceRand<double>(0, 1).Entropy()), WithinRel(1 + std::log(2.0), 1e-15));
    CHECK_THAT(static_cast<double>(LaplaceRand<double>(4, 2).Entropy()), WithinRel(1 + std::log(4.0), 1e-15));
    CHECK_THAT(static_cast<double>(AsymmetricLaplaceRand<double>(7, 3, 2).Entropy()),
               WithinRel(1 + std::log(7.5), 1e-15));
}

TEST_CASE("zero noncentrality is the central chi-squared limit", "[regression]")
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

    std::mt19937_64 actualGenerator(42), expectedGenerator(42);
    const double actual = NoncentralChiSquaredRand<double>::Variate(0.5, 0, actualGenerator);
    const double expected = 2 * GammaDistribution<double>::StandardVariate(0.25, expectedGenerator);
    CHECK(actual == expected);
}

TEST_CASE("Bernoulli endpoint entropy uses the zero-log-zero limit", "[regression]")
{
    CHECK(BernoulliRand(0).Entropy() == 0);
    CHECK(BernoulliRand(1).Entropy() == 0);
    CHECK_THAT(BernoulliRand(0.5).Entropy(), WithinRel(std::log(2.0), 1e-15));
}

TEST_CASE("binomial probability boundaries preserve construction and sampling", "[regression]")
{
    const BinomialRand<int> certain(7, 1);
    std::mt19937_64 certainGenerator(42);
    CHECK(certain.Mode() == 7);
    CHECK(certain(certainGenerator) == 7);
    std::vector<int> sample(3);
    certain.Sample(certainGenerator, sample);
    CHECK(sample == std::vector<int>{7, 7, 7});

    const BinomialRand<int> impossible(7, 0);
    std::mt19937_64 impossibleGenerator(42);
    CHECK(impossible.Mode() == 0);
    CHECK(impossible(impossibleGenerator) == 0);
    CHECK(BernoulliRand(1).Mode() == 1);
    CHECK(BinomialRand<int>(5, 0.5).Mode() == 3);

    BinomialRand<int> bulk(10000000, 1 - 9e-7), scalar(10000000, 1 - 9e-7);
    std::mt19937_64 bulkGenerator(42), scalarGenerator(42);
    std::vector<int> actual(4), expected(4);
    bulk.Sample(bulkGenerator, actual);
    for (int &value : expected)
        value = scalar(scalarGenerator);
    CHECK(actual == expected);
}

TEST_CASE("discrete-uniform characteristic function handles lattice periods", "[regression]")
{
    const UniformDiscreteRand<int> distribution(-3, 7);
    const double period = 2 * static_cast<double>(randlib::consts::pi);
    CHECK_THAT(distribution.CF(period).real(), WithinAbs(1.0, 0.0));
    CHECK_THAT(distribution.CF(period).imag(), WithinAbs(0.0, 0.0));
    CHECK_THAT(distribution.CF(2 * period).real(), WithinAbs(1.0, 0.0));
    CHECK_THAT(distribution.CF(-period).imag(), WithinAbs(0.0, 0.0));
}

TEST_CASE("Skellam mode maximizes the PMF instead of truncating the mean", "[regression]")
{
    const SkellamRand<int> distribution(2, 1.1);
    CHECK_THAT(distribution.P(0), WithinRel(0.2139951994197166138511878227306453, 1e-13));
    CHECK_THAT(distribution.P(1), WithinRel(0.2329981383802911800024090128922231, 1e-13));
    CHECK_THAT(distribution.P(2), WithinRel(0.1772656913264927706363333023355158, 1e-13));
    CHECK(distribution.Mode() == 1);
}

TEST_CASE("static Poisson sampling uses the exact cached algorithm", "[regression]")
{
    for (double rate : {5.0, 1001.0}) {
        const PoissonRand<int> distribution(rate);
        std::mt19937_64 staticGenerator(42), instanceGenerator(42);
        CHECK(PoissonRand<int>::Variate(rate, staticGenerator) == distribution(instanceGenerator));
    }
}

TEST_CASE("negative-binomial success endpoint is deterministic", "[regression]")
{
    PascalRand<int> negativeBinomial(3, 1);
    CHECK(negativeBinomial.P(0) == 1);
    CHECK(negativeBinomial.logP(0) == 0);
    CHECK(negativeBinomial.P(1) == 0);
    CHECK(negativeBinomial.logP(1) == -INFINITY);
    CHECK(negativeBinomial.F(0) == 1);
    CHECK(negativeBinomial.S(0) == 0);
    CHECK(negativeBinomial.Mean() == 0);
    CHECK(negativeBinomial.Variance() == 0);
    CHECK(std::isnan(negativeBinomial.Skewness()));
    CHECK(std::isnan(negativeBinomial.ExcessKurtosis()));
    CHECK(negativeBinomial.Mode() == 0);
    std::mt19937_64 negativeBinomialGenerator(42);
    CHECK(negativeBinomial(negativeBinomialGenerator) == 0);
    std::vector<int> sample(3, -1);
    negativeBinomial.Sample(negativeBinomialGenerator, sample);
    CHECK(sample == std::vector<int>{0, 0, 0});

    GeometricRand<int> geometric(1);
    CHECK(geometric.P(0) == 1);
    CHECK(geometric.logP(0) == 0);
    CHECK(geometric.P(1) == 0);
    CHECK(geometric.logP(1) == -INFINITY);
    CHECK(geometric.F(0) == 1);
    CHECK(geometric.S(0) == 0);
    CHECK(geometric.Median() == 0);
    std::mt19937_64 generator(42);
    CHECK(geometric(generator) == 0);
    CHECK(GeometricRand<int>::Variate(1, generator) == 0);
    CHECK_THROWS_AS(GeometricRand<int>::Variate(0, generator), std::invalid_argument);
}

TEST_CASE("geometric entropy uses natural units", "[regression]")
{
    CHECK(GeometricRand<int>(1).Entropy() == 0);
    CHECK_THAT(static_cast<double>(GeometricRand<int>(0.5).Entropy()), WithinRel(2 * std::log(2.0), 1e-15));
}

TEST_CASE("bivariate covariance uses marginal standard deviations", "[regression]")
{
    const BivariateNormalRand<double> normal(1, -2, 2, 3, 0.25);
    const LongDoubleTriplet normalCovariance = normal.Covariance();
    CHECK_THAT(static_cast<double>(std::get<0>(normalCovariance)), WithinRel(4.0, 1e-15));
    CHECK_THAT(static_cast<double>(std::get<1>(normalCovariance)), WithinRel(1.5, 1e-15));
    CHECK_THAT(static_cast<double>(std::get<2>(normalCovariance)), WithinRel(9.0, 1e-15));

    const TrinomialRand<int> trinomial(10, 0.2, 0.3);
    CHECK_THAT(static_cast<double>(trinomial.Correlation()), WithinRel(-std::sqrt(3.0 / 28.0), 1e-15));
    const LongDoubleTriplet trinomialCovariance = trinomial.Covariance();
    CHECK_THAT(static_cast<double>(std::get<1>(trinomialCovariance)), WithinAbs(-0.6, 1e-15));

    const TrinomialRand<int> degenerate(5, 0, 0.4);
    CHECK(std::isnan(degenerate.Correlation()));
    CHECK(std::get<1>(degenerate.Covariance()) == 0);
}

TEST_CASE("normal-inverse-gamma marginals, sampling and joint CDF agree", "[regression]")
{
    NormalInverseGammaRand<double> distribution(2, 3, 4, 2);
    const auto marginals = distribution.GetMarginalDistributions();
    CHECK_THAT(marginals.first.GetScale(), WithinRel(std::sqrt(1.0 / 6.0), 1e-15));
    const LongDoubleTriplet covariance = distribution.Covariance();
    CHECK_THAT(static_cast<double>(std::get<0>(covariance)), WithinRel(2.0 / 9.0, 1e-15));
    CHECK(std::get<1>(covariance) == 0);
    CHECK_THAT(static_cast<double>(std::get<2>(covariance)), WithinRel(2.0 / 9.0, 1e-15));

    const double sigmaSq = 1.0;
    CHECK_THAT(distribution.F({2.0, sigmaSq}), WithinRel(0.5 * marginals.second.F(sigmaSq), 1e-15));
    CHECK_THAT(distribution.F({1.0, sigmaSq}), WithinAbs(0.009262835921861568182061857402, 5e-11));
    CHECK_THAT(distribution.F({3.0, sigmaSq}), WithinAbs(0.8478606245766854804799349441, 5e-11));
    CHECK(distribution.F({-INFINITY, sigmaSq}) == 0);
    CHECK(distribution.F({INFINITY, sigmaSq}) == marginals.second.F(sigmaSq));
    CHECK(distribution.F({1.0, INFINITY}) == marginals.first.F(1.0));

    std::mt19937_64 generator(42);
    std::mt19937_64 expectedGenerator(42);
    InverseGammaRand<double> inverseGamma(4, 2);
    const double expectedSigmaSq = inverseGamma(expectedGenerator);
    const double expectedX = 2 + std::sqrt(expectedSigmaSq / 3) * detail::standardNormal<double>(expectedGenerator);
    CHECK(distribution(generator) == std::make_pair(expectedX, expectedSigmaSq));

    CHECK(std::isnan(NormalInverseGammaRand<double>(0, 1, 1, 1).Correlation()));
}

TEST_CASE("negative-correlation bivariate normal CDF uses standardized limits", "[regression]")
{
    const BivariateNormalRand<double> distribution(1, -2, 2, 0.5, -1);
    CHECK_THAT(distribution.F({3.0, -2.0}), WithinAbs(0.3413447460685429485852325456, 1e-15));
    CHECK(distribution.F({1.0, -3.0}) == 0);
}

TEST_CASE("trinomial endpoints and modes stay on the probability simplex", "[regression]")
{
    const TrinomialRand<int> boundary(5, 0.4, 0.6);
    CHECK_THAT(boundary.P({2, 3}), WithinAbs(0.3456, 1e-15));
    CHECK_THAT(boundary.F({2, 3}), WithinAbs(0.3456, 1e-15));
    CHECK(std::isfinite(boundary.logP({2, 3})));
    CHECK(boundary.logP({1, 1}) == -INFINITY);

    TrinomialRand<int> firstDeterministic(5, 1, 0);
    CHECK(firstDeterministic.P({5, 0}) == 1);
    CHECK(firstDeterministic.logP({5, 0}) == 0);
    CHECK(firstDeterministic.logP({0, 0}) == -INFINITY);
    CHECK(firstDeterministic.F({4, 5}) == 0);
    CHECK(firstDeterministic.F({5, 0}) == 1);
    std::mt19937_64 firstGenerator(42);
    CHECK(firstDeterministic(firstGenerator) == std::make_pair(5, 0));

    TrinomialRand<int> secondDeterministic(5, 0, 1);
    CHECK(secondDeterministic.P({0, 5}) == 1);
    CHECK(secondDeterministic.logP({0, 5}) == 0);
    CHECK(secondDeterministic.logP({0, 0}) == -INFINITY);
    CHECK(secondDeterministic.F({5, 4}) == 0);
    CHECK(secondDeterministic.F({0, 5}) == 1);
    std::mt19937_64 secondGenerator(42);
    CHECK(secondDeterministic(secondGenerator) == std::make_pair(0, 5));

    const TrinomialRand<int> formerInvalidMode(1, 0.26, 0.49);
    CHECK(formerInvalidMode.Mode() == std::make_pair(0, 1));
    CHECK(formerInvalidMode.Mode().first + formerInvalidMode.Mode().second <= 1);

    const std::vector<double> probabilities{0.0, 0.1, 0.25, 0.4, 0.6, 1.0};
    for (int n = 1; n <= 6; ++n) {
        for (double p1 : probabilities) {
            for (double p2 : probabilities) {
                if (p1 + p2 > 1.0)
                    continue;
                const TrinomialRand<int> distribution(n, p1, p2);
                const auto mode = distribution.Mode();
                double maximum = 0.0;
                for (int x = 0; x <= n; ++x)
                    for (int y = 0; y <= n - x; ++y)
                        maximum = std::max(maximum, distribution.P({x, y}));
                CAPTURE(n, p1, p2, mode);
                CHECK(mode.first >= 0);
                CHECK(mode.second >= 0);
                CHECK(mode.first + mode.second <= n);
                CHECK_THAT(distribution.P(mode), WithinAbs(maximum, 1e-15));
            }
        }
    }
}

TEST_CASE("hypergeometric moments cover small, deterministic and large populations", "[regression]")
{
    const HyperGeometricRand<int> singleton(1, 1, 1);
    CHECK(singleton.Mean() == 1);
    CHECK(singleton.Variance() == 0);
    CHECK(std::isnan(singleton.Skewness()));
    CHECK(std::isnan(singleton.ExcessKurtosis()));

    const HyperGeometricRand<int> symmetricBernoulli(2, 1, 1);
    CHECK(symmetricBernoulli.Mean() == 0.5);
    CHECK(symmetricBernoulli.Variance() == 0.25);
    CHECK(symmetricBernoulli.Skewness() == 0);
    CHECK_THAT(static_cast<double>(symmetricBernoulli.ExcessKurtosis()), WithinAbs(-2.0, 1e-15));

    const HyperGeometricRand<int> smallAsymmetric(3, 1, 1);
    CHECK_THAT(static_cast<double>(smallAsymmetric.Skewness()), WithinAbs(1.0 / std::sqrt(2.0), 1e-15));
    CHECK_THAT(static_cast<double>(smallAsymmetric.ExcessKurtosis()), WithinAbs(-1.5, 1e-15));

    const HyperGeometricRand<int> deterministic(10, 10, 4);
    CHECK(deterministic.Mean() == 4);
    CHECK(deterministic.Variance() == 0);
    CHECK(std::isnan(deterministic.Skewness()));
    CHECK(std::isnan(deterministic.ExcessKurtosis()));

    const HyperGeometricRand<int> large(2000000000, 1500000000, 1000000000);
    CHECK(large.Mean() == 750000000);
    CHECK(large.Mode() == 750000000);

    const HyperGeometricRand<int> shiftedSupport(100, 80, 70);
    CHECK_THAT(shiftedSupport.F(shiftedSupport.MinValue()),
               WithinAbs(shiftedSupport.P(shiftedSupport.MinValue()), 1e-15));
}

TEST_CASE("negative-hypergeometric moments do not overflow before conversion", "[regression]")
{
    const int population = std::numeric_limits<int>::max();
    const int successes = population - 1;
    const int target = 1000000000;
    const NegativeHyperGeometricRand<int> distribution(population, successes, target);
    const long double expectedMean = static_cast<long double>(target) / population;
    const long double expectedVariance = target * (1.0L - static_cast<long double>(target) / population) / population;
    CHECK_THAT(static_cast<double>(distribution.Mean()), WithinRel(static_cast<double>(expectedMean), 1e-15));
    CHECK_THAT(static_cast<double>(distribution.Variance()), WithinRel(static_cast<double>(expectedVariance), 1e-15));
    CHECK(distribution.Variance() > 0);
}

TEST_CASE("beta-family distribution functions are complementary at the ends of the support", "[regression]")
{
    /// 1 - x rounds to 1 for x < 1e-16, which made S(x) = 1 while F(x) = 0.29
    const BetaRand<double> beta(0.0175, 0.442);
    CHECK_THAT(beta.F(1e-30), WithinRel(0.28995773045471592299, 1e-14));
    CHECK_THAT(beta.S(1e-30), WithinRel(0.71004226954528407701, 1e-14));
    /// x / (1 + x) rounds to 1 for large x, which made F(x) = 1 while S(x) = 1e-4
    const BetaPrimeRand<double> betaPrime(8.97, 0.234);
    CHECK_THAT(betaPrime.F(1.5e18), WithinRel(0.99989851563113950725, 1e-14));
    CHECK_THAT(betaPrime.S(1.5e18), WithinRel(0.00010148436886049275311, 1e-13));
    CHECK(betaPrime.F(INFINITY) == 1);
    CHECK(betaPrime.S(INFINITY) == 0);
    for (double shape : {0.281442, 0.5}) {
        const ArcsineRand<double> arcsine(shape, -2, 5);
        for (double x : {-1.9999999999, 4.9999999, std::nextafter(5.0, 0.0)})
            CHECK_THAT(arcsine.F(x) + arcsine.S(x), WithinAbs(1.0, 2e-16));
    }
}

TEST_CASE("Irwin-Hall upper tail uses the symmetric lower-half sum", "[regression]")
{
    const IrwinHallRand<double> distribution(28);
    CHECK_THAT(distribution.S(26), WithinRel(8.804384711454354053e-22, 1e-12));
    CHECK_THAT(distribution.S(24.5), WithinRel(5.6076188672235636108e-15, 1e-12));
    CHECK(distribution.f(28) == 0);
    CHECK_THAT(distribution.f(20), WithinRel(distribution.f(8), 1e-15));
}

TEST_CASE("Zeta survival function does not subtract from one", "[regression]")
{
    const ZetaRand<int> distribution(5);
    CHECK_THAT(distribution.S(20), WithinRel(1.3624373201755622612e-6, 1e-14));
    CHECK_THAT(distribution.S(1000), WithinRel(2.4061504326487821659e-13, 1e-14));
    CHECK_THAT(distribution.F(15) + distribution.S(15), WithinAbs(1.0, 5e-16));
}

TEST_CASE("Planck survival function integrates the upper tail", "[regression]")
{
    const PlanckRand<double> distribution(12.8609, 1.45836);
    CHECK_THAT(distribution.S(35), WithinRel(1.9193573707032127246e-10, 1e-12));
}

TEST_CASE("Skellam probabilities beyond the range of the Bessel function", "[regression]")
{
    const SkellamRand<int> distribution(190.832, 0.203018);
    CHECK_THAT(distribution.P(275), WithinRel(1.8628462183563493992e-9, 1e-12));
    CHECK_THAT(distribution.P(290), WithinRel(5.0207397258500101208e-12, 1e-12));
}

TEST_CASE("exponentially modified Gaussian lower tail with a large rate", "[regression]")
{
    /// erfc underflowed while the exponential factor overflowed: f and F returned NaN
    const ExponentiallyModifiedGaussianRand<double> distribution(-7.65171, 62.6909, 10.4122);
    CHECK_THAT(distribution.f(-50), WithinRel(2.9046476425714958643e-8, 1e-13));
    CHECK_THAT(distribution.F(-50), WithinRel(4.1548525164977641925e-8, 1e-13));
    CHECK_THAT(distribution.Quantile(1e-10), WithinRel(-57.926758974639725, 1e-12));
}

TEST_CASE("generalized extreme-value moments near the Gumbel limit", "[regression]")
{
    const GeneralizedExtremeValueRand<double> nearGumbel(1e-8, 2, 1);
    CHECK_THAT(static_cast<double>(nearGumbel.Mean()), WithinRel(2.1544313495841858093, 1e-15));
    CHECK_THAT(static_cast<double>(nearGumbel.Variance()), WithinRel(6.5797364395159995752, 1e-14));
    CHECK_THAT(static_cast<double>(nearGumbel.Skewness()), WithinRel(1.1395471590707740121, 1e-13));
    CHECK_THAT(static_cast<double>(nearGumbel.ExcessKurtosis()), WithinRel(2.4000002891768986698, 1e-13));
    const GeneralizedExtremeValueRand<double> bounded(-1e-4, 2, 1);
    CHECK_THAT(static_cast<double>(bounded.Skewness()), WithinRel(1.1389505609250348749, 1e-13));
    CHECK_THAT(static_cast<double>(bounded.ExcessKurtosis()), WithinRel(2.3971097566600894297, 1e-13));
}

TEST_CASE("negative-binomial sampling with a small shape", "[regression]")
{
    /// the mixing gamma variate underflows to zero, the limit of a zero Poisson variate
    NegativeBinomialRand<int, double> distribution(0.02, 0.2);
    distribution.Reseed(3);
    std::vector<int> sample(200000);
    REQUIRE_NOTHROW(distribution.Sample(sample));
    CHECK(std::ranges::all_of(sample, [](int value) { return value >= 0; }));
}

TEST_CASE("characteristic functions of narrow laws far from the origin", "[regression]")
{
    auto distance = [](std::complex<double> actual, double re, double im) {
        return std::abs(actual - std::complex<double>(re, im));
    };
    CHECK(distance(FrechetRand<double>(18.4879, 2.91585, -1.59155).CF(5.093114), 0.45384760893887839673,
                   0.43044296347465515653) < 1e-12);
    CHECK(distance(RiceRand<double>(3.53436, 0.11411).CF(8.357852), -0.18155087820694433489, -0.60821636381094294381) <
          1e-12);
    CHECK(distance(GumbelRand<double>(100, 0.01).CF(30), -0.92900728149985855152, 0.054095335821230002033) < 1e-12);
}

TEST_CASE("Rice moments for a large offset-to-scale ratio", "[regression]")
{
    /// the Bessel functions in the mean overflow; their Hankel expansion keeps the mean exact
    const RiceRand<double> distribution(50, 0.1);
    CHECK_THAT(static_cast<double>(distribution.Mean()), WithinRel(50.0001000001000006, 1e-15));
    CHECK_THAT(static_cast<double>(distribution.Variance()), WithinRel(0.00999997999991999912, 1e-9));
    CHECK_THAT(static_cast<double>(distribution.Skewness()), WithinRel(8.0000960016200387772e-9, 1e-3));
}

TEST_CASE("log-normal characteristic function keeps its relative accuracy", "[regression]")
{
    auto relative = [](std::complex<double> actual, double re, double im) {
        const std::complex<double> expected(re, im);
        return std::abs(actual - expected) / std::abs(expected);
    };
    CHECK(relative(LogNormalRand<double>(0, 0.01).CF(100), -3.9894565854085051679e-13, -5.0063958950743116855e-13) <
          1e-12);
    CHECK(relative(LogNormalRand<double>(3, 0.01).CF(10), 8.2913213202566704761e-29, 3.0916895064135382773e-31) <
          1e-12);
    CHECK(relative(LogNormalRand<double>(0, 16).CF(1), 0.44363937308484652105, 0.15189506168865027344) < 1e-12);
}

TEST_CASE("bulk evaluations fill no more values than the output holds", "[regression]")
{
    const NormalRand<double> normal;
    const std::vector<double> x = {-1, 0, 1};
    std::vector<double> y(2, -1);
    normal.ProbabilityDensityFunction(x, y);
    CHECK(y == std::vector<double>{normal.f(-1), normal.f(0)});
    normal.CumulativeDistributionFunction(x, y);
    CHECK(y == std::vector<double>{normal.F(-1), normal.F(0)});
    std::vector<double> quantiles(1);
    normal.QuantileFunction(std::vector<double>{0.5, 0.9}, quantiles);
    CHECK(quantiles[0] == 0);

    const PoissonRand<int> poisson(2);
    poisson.LogProbabilityMassFunction(std::vector<int>{0, 1, 2}, y);
    CHECK(y == std::vector<double>{poisson.logP(0), poisson.logP(1)});
}

TEST_CASE("incomplete beta and stable quadrature avoid their slowest regimes", "[regression]")
{
    /// the incomplete beta series converge at the rate x near 1; maximum likelihood probes such degrees
    const auto start = std::chrono::steady_clock::now();
    const NoncentralTRand<double> tinyDegree(2.4598888848425963e-07, 40.274620841745602);
    const NoncentralTRand<double> smallDegree(2.7885942799549225e-06, 26.469118041391997);
    CHECK(std::isfinite(tinyDegree.logf(11)));
    CHECK(std::isfinite(smallDegree.logf(1.659707541784905)));
    /// nested stable integrals whose adaptive rule cannot reach the tolerance stop at the rounding floor
    const StableCountRand<double> stableCount(0.498942);
    CHECK_THAT(stableCount.S(20) + stableCount.F(20), WithinAbs(1.0, 1e-12));
    CHECK(stableCount.S(50) > 0);
    CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(2));
}

TEST_CASE("von Mises tails keep their relative accuracy", "[regression]")
{
    /// the Fourier series has an absolute error of a few ulps, which used to dominate both small tails
    const VonMisesRand<double> distribution(0, 6.59014);
    CHECK_THAT(distribution.F(-3.1), WithinRel(7.8900301649572674624e-8, 1e-12));
    CHECK_THAT(distribution.S(2.9), WithinRel(4.8842392821740510472e-7, 1e-12));
    CHECK(distribution.F(-1.0) + distribution.S(-1.0) == 1.0);
}

TEST_CASE("Brent minimization terminates far from the origin", "[regression]")
{
    /// an absolute tolerance below the spacing of doubles near 1e9 never ended the iteration
    double root = 0;
    math::findMin<double>([](double x) { return (x - 1e9) * (x - 1e9); }, 1e9 + 5, root);
    CHECK_THAT(root, WithinRel(1e9, 1e-12));

    /// a mean of 1e10 used to start the search for a mode of order one
    const StableCountRand<double> stableCount(0.110525, 0.603312, -1.35676);
    const double mode = stableCount.Mode();
    const double step = 1e-4;
    CHECK(stableCount.logf(mode) >= stableCount.logf(mode - step));
    CHECK(stableCount.logf(mode) >= stableCount.logf(mode + step));
    CHECK_THAT(StableCountRand<double>(0.5, 2, 1).Mode(), WithinRel(5.0, 1e-7));
}
