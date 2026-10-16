#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/continuous/log_stable.hpp>
#include <randlib/continuous/stable.hpp>

#include <array>
#include <cmath>
#include <span>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using randlib::LogStableRand;
using randlib::Pcg64;
using randlib::StableRand;

TEST_CASE("stable S0 and S1 parameterizations describe the same law", "[stable]")
{
    using Parameterization = StableRand<double>::Parameterization;
    constexpr double pi = 3.14159265358979323846;

    for (const auto &[alpha, beta, scale, location0] :
         std::array<std::array<double, 4>, 3>{{{0.7, -0.4, 1.5, 2.0}, {1.0, 0.6, 2.0, -1.0}, {1.4, 0.3, 0.75, 4.0}}}) {
        const double shift =
            alpha == 1.0 ? 2.0 / pi * beta * scale * std::log(scale) : beta * scale * std::tan(0.5 * pi * alpha);
        StableRand<double> s0(alpha, beta, scale, location0, Parameterization::S0);
        StableRand<double> s1(alpha, beta, scale, location0 - shift, Parameterization::S1);

        CHECK(s0.GetParameterization() == Parameterization::S0);
        CHECK(s0.GetLocation() == location0);
        CHECK(s1.GetLocation() == location0 - shift);
        for (double x : {-3.0, 0.0, 5.0}) {
            CHECK_THAT(s0.f(x), WithinRel(s1.f(x), 2e-14));
            CHECK_THAT(s0.F(x), WithinAbs(s1.F(x), 2e-14));
            CHECK_THAT(s0.S(x), WithinAbs(s1.S(x), 2e-14));
        }
        for (double t : {-2.0, 0.5, 3.0})
            CHECK_THAT(std::abs(s0.CF(t) - s1.CF(t)), WithinAbs(0.0, 2e-14));

        Pcg64 first(42), second(42);
        for (int i = 0; i < 32; ++i)
            CHECK_THAT(s0(first), WithinRel(s1(second), 2e-14));
    }
}

TEST_CASE("stable S0 location has the expected support and setter semantics", "[stable]")
{
    using Parameterization = StableRand<double>::Parameterization;
    StableRand<double> stable(0.5, 1.0, 2.0, 3.0, Parameterization::S0);

    CHECK(stable.GetLocation() == 3.0);
    CHECK_THAT(stable.MinValue(), WithinAbs(1.0, 1e-14));
    CHECK(stable.Name() == "Stable(0.5, 1, 2, 3, S0)");

    stable.SetScale(4.0);
    CHECK(stable.GetLocation() == 3.0);
    CHECK_THAT(stable.MinValue(), WithinAbs(-1.0, 1e-14));
    stable.SetLocation(7.0);
    CHECK(stable.GetLocation() == 7.0);
    CHECK_THAT(stable.MinValue(), WithinAbs(3.0, 1e-14));

    const auto characteristicFunction = stable.CF(0.75);
    stable.SetParameterization(Parameterization::S1);
    CHECK_THAT(stable.GetLocation(), WithinAbs(3.0, 1e-14));
    CHECK_THAT(std::abs(stable.CF(0.75) - characteristicFunction), WithinAbs(0.0, 1e-14));
    CHECK(stable.Name() == "Stable(0.5, 1, 4, 3)");
}

TEST_CASE("log-stable probability and sampling are exponential stable transforms", "[stable][log-stable]")
{
    using LogStable = LogStableRand<double>;
    const LogStable distribution(1.5, -1, 0.5, 0.25);
    const StableRand<double> stable(1.5, -1, 0.5, 0.25);
    for (double x : {0.1, 0.5, 1.0, 2.0, 10.0}) {
        CHECK_THAT(distribution.f(x), WithinRel(stable.f(std::log(x)) / x, 2e-14));
        CHECK_THAT(distribution.logf(x), WithinAbs(stable.logf(std::log(x)) - std::log(x), 2e-14));
        CHECK_THAT(distribution.F(x), WithinAbs(stable.F(std::log(x)), 2e-14));
        CHECK_THAT(distribution.S(x), WithinAbs(stable.S(std::log(x)), 2e-14));
    }
    for (double probability : {0.01, 0.3, 0.7, 0.99}) {
        CHECK_THAT(distribution.Quantile(probability), WithinRel(std::exp(stable.Quantile(probability)), 2e-14));
        CHECK_THAT(distribution.Quantile1m(probability), WithinRel(std::exp(stable.Quantile1m(probability)), 2e-14));
    }

    Pcg64 first(42), second(42);
    for (int i = 0; i < 32; ++i)
        CHECK(distribution(first) == std::exp(stable(second)));

    LogStable s0 = distribution;
    s0.SetParameterization(LogStable::Parameterization::S0);
    CHECK(s0.GetParameterization() == LogStable::Parameterization::S0);
    CHECK_THAT(s0.F(1), WithinAbs(distribution.F(1), 0));
    CHECK_THAT(s0.Mean(), WithinAbs(distribution.Mean(), 0));
}

TEST_CASE("log-stable moments follow the spectrally negative Laplace exponent", "[stable][log-stable]")
{
    struct Reference
    {
        LogStableRand<double> distribution;
        double mean;
        double variance;
        double skewness;
        double excess;
    };
    const std::array references{
        Reference{LogStableRand<double>(0.5, -1, 1, 0), 0.24311673443421421, 0.076229536674656454, 0.77779277586209408,
                  -0.82146413679350378},
        Reference{LogStableRand<double>(1, -1, 0.5, 0.25), 1.2840254166877415, 0.9145164660676866, 0.46240949219948345,
                  -0.56276865611099966},
        Reference{LogStableRand<double>(1.5, -1, 0.5, 0.25), 2.1170000166126747, 2.2999143208767069, 1.2481512140756379,
                  2.5991895948497396},
        Reference{LogStableRand<double>(2, 0, 0.5, 0), 1.2840254166877415, 1.0695605577589171, 2.9387981481747263,
                  18.507279724983916},
    };
    for (const auto &[distribution, mean, variance, skewness, excess] : references) {
        CHECK_THAT(distribution.Mean(), WithinRel(mean, 2e-14));
        CHECK_THAT(distribution.Variance(), WithinRel(variance, 2e-14));
        CHECK_THAT(distribution.Skewness(), WithinRel(skewness, 2e-13));
        CHECK_THAT(distribution.ExcessKurtosis(), WithinRel(excess, 2e-13));
    }

    const LogStableRand<double> withoutMoments(1.5, 0, 1, 0);
    CHECK(std::isinf(withoutMoments.Mean()));
    CHECK(std::isinf(withoutMoments.Variance()));
    CHECK(std::isnan(withoutMoments.Skewness()));
    CHECK(std::isnan(withoutMoments.ExcessKurtosis()));

    const LogStableRand<double> concentratedGaussian(2, 0, 1e-10, 0);
    CHECK_THAT(concentratedGaussian.Skewness(), WithinRel(3 * std::sqrt(2.0) * 1e-10, 2e-14));
    CHECK_THAT(concentratedGaussian.ExcessKurtosis(), WithinRel(32e-20, 2e-14));
}

TEST_CASE("log-stable mode accounts for the exponential Jacobian", "[stable][log-stable]")
{
    const LogStableRand<double> gaussian(2, 0, 0.5, 0);
    CHECK_THAT(gaussian.Median(), WithinAbs(1, 1e-14));
    CHECK_THAT(gaussian.Mode(), WithinRel(std::exp(-0.5), 1e-14));

    const LogStableRand<double> heavyLeftTail(1.5, -1, 0.5, 0.25);
    CHECK(heavyLeftTail.Mode() == 0);
    CHECK(std::isinf(heavyLeftTail.f(0)));

    const LogStableRand<double> levy(0.5, 1, 1, 0);
    const double logMode = 0.25 * (std::sqrt(17.0) - 3);
    CHECK_THAT(levy.Mode(), WithinRel(std::exp(logMode), 2e-8));
    CHECK(levy.MinValue() == 1);

    CHECK(gaussian.CF(0) == std::complex<double>(1, 0));
    CHECK_THROWS_AS(gaussian.CF(1), randlib::not_implemented);
}

TEST_CASE("log-stable short-tail quantiles preserve the finite skew-angle endpoint", "[stable][log-stable]")
{
    const LogStableRand<double> distribution(1.5, -1, 0.5, 0.25);
    CHECK_THAT(distribution.Quantile1m(1e-10), WithinRel(33.17281190476198487520755, 1e-10));
}

TEST_CASE("stable probability functions cover the former alpha-one gap", "[stable]")
{
    // Generated by tools/refdata/distributions.py with mpmath at 40 decimal digits.
    struct Reference
    {
        double x;
        double density;
        double cdf;
    };
    const std::array<Reference, 5> nearOne{{
        {-5, 0.0000333588035793094259245322585947, 0.00233254382039488143159070939797},
        {-1, 0.0000375528674040140126750071980457, 0.00247411828623862250083135373738},
        {0, 0.0000387253487862958344983538734164, 0.00251225287729380940184560462047},
        {1, 0.0000399537567485926291256994478015, 0.00255158762299297559911006656114},
        {5, 0.0000455037248943888733862200239612, 0.00272214083768517602308541890061},
    }};
    const std::array<Reference, 5> smallSkew{{
        {-5, 0.0121752297462174968836042920268, 0.062479570601748421299710799891},
        {-1, 0.159194375963569286555405426159, 0.249134143625264197000340978652},
        {0, 0.318305866405370471170519969017, 0.499415152353592157861622646126},
        {1, 0.159117147383519696630710656334, 0.749134088960657294736048344699},
        {5, 0.012310197475911957547551468963, 0.936813264635769399039323483972},
    }};

    const auto check = [](const StableRand<double> &distribution, std::span<const Reference> reference) {
        for (const auto &[x, density, cdf] : reference) {
            CHECK_THAT(distribution.f(x), WithinAbs(density, 1e-9));
            CHECK_THAT(distribution.F(x), WithinAbs(cdf, 1e-9));
            CHECK_THAT(distribution.S(x), WithinAbs(1 - cdf, 1e-9));
        }
    };
    check(StableRand<double>(0.995, 0.5, 1, 0), nearOne);
    check(StableRand<double>(1, 0.005, 1, 0), smallSkew);
}

TEST_CASE("stable S0 probability functions are continuous across alpha one", "[stable]")
{
    // Generated by tools/refdata/distributions.py with mpmath at 40 decimal digits.
    using Parameterization = StableRand<double>::Parameterization;
    struct Reference
    {
        double alpha;
        double x;
        double density;
        double cdf;
    };
    const std::array<Reference, 15> references{{
        {0.995, -5, 0.005797400803804647049855652, 0.02972034690482593868278839},
        {0.995, 0, 0.2926530326779768599443813, 0.4373257444117281642898867},
        {0.995, 5, 0.01927169947761192167594759, 0.8989530504911191519663653},
        {0.999, -5, 0.005782469239752549147422102, 0.02949255366561870614446809},
        {0.999, 0, 0.2925468938522064997428655, 0.4374742218884701834735273},
        {0.999, 5, 0.01923155890825204956342629, 0.8996930075144387738015242},
        {1.0, -5, 0.005778717168275644635642541, 0.02943582495669952017340078},
        {1.0, 0, 0.2925204705660767133376541, 0.4375114838590878790180263},
        {1.0, 5, 0.01922144475475842389376862, 0.8998773763860317059289712},
        {1.001, -5, 0.005774957470151258645376014, 0.02937918380520343741977727},
        {1.001, 0, 0.2924940918270732129513717, 0.4375488026945080022459222},
        {1.001, 5, 0.01921129912827762126682471, 0.9000614977613298338018618},
        {1.005, -5, 0.005759842863944214185283139, 0.02915349195569516979912785},
        {1.005, 0, 0.292389022124684687472877, 0.4376986449438168502457942},
        {1.005, 5, 0.01917040320174463719899932, 0.9007955148513268842718892},
    }};

    for (const auto &[alpha, x, density, cdf] : references) {
        const StableRand<double> distribution(alpha, 0.5, 1, 0, Parameterization::S0);
        CHECK_THAT(distribution.f(x), WithinAbs(density, 1e-9));
        CHECK_THAT(distribution.F(x), WithinAbs(cdf, 1e-9));
        CHECK_THAT(distribution.S(x), WithinAbs(1 - cdf, 1e-9));
    }

    const StableRand<double> atOne(1, 0.5, 1, 0, Parameterization::S0);
    for (double alpha :
         {1 - 1e-8, 1 + 1e-8, 1 - 1e-12, 1 + 1e-12, std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0)}) {
        const StableRand<double> adjacent(alpha, 0.5, 1, 0, Parameterization::S0);
        CHECK_THAT(adjacent.f(0), WithinAbs(atOne.f(0), 1e-8));
        CHECK_THAT(adjacent.F(0), WithinAbs(atOne.F(0), 1e-8));
        CHECK_THAT(std::abs(adjacent.CF(0.5) - atOne.CF(0.5)), WithinAbs(0.0, 1e-8));
    }

    for (double alpha : {1 - 1e-12, 1 + 1e-12, std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0)}) {
        Pcg64 exactEngine(19), adjacentEngine(19);
        const StableRand<double> adjacent(alpha, 0.5, 1, 0, Parameterization::S0);
        for (int i = 0; i < 16; ++i)
            CHECK_THAT(adjacent(adjacentEngine), WithinAbs(atOne(exactEngine), 1e-6));
    }

    Pcg64 exactBulkEngine(23), adjacentBulkEngine(23);
    std::array<double, 16> exactBulk{}, adjacentBulk{};
    const StableRand<double> adjacent(std::nextafter(1.0, 2.0), 0.5, 1, 0, Parameterization::S0);
    atOne.Sample(exactBulkEngine, exactBulk);
    adjacent.Sample(adjacentBulkEngine, adjacentBulk);
    for (std::size_t i = 0; i < exactBulk.size(); ++i)
        CHECK_THAT(adjacentBulk[i], WithinAbs(exactBulk[i], 1e-6));

    for (double probability : {0.01, 0.3, 0.7, 0.99}) {
        const double quantile = adjacent.Quantile(probability);
        CHECK_THAT(adjacent.F(quantile), WithinAbs(probability, 2e-10));
    }

    struct TailReference
    {
        double x;
        double density;
        double cdf;
    };
    const std::array<TailReference, 5> farReferences{{
        {100, 0.00004887150862277182256299117, 0.9951614753933674957184178},
        {1000, 0.000000479287647385920136523552, 0.9995215481785292461544734},
        {10000, 0.000000004777168222832025456243, 0.9999522401586044993935791},
        {100000, 0.00000000004774970208682687927025, 0.9999952251831522225273105},
        {1000000, 0.0000000000004774687482163897010552, 0.9999995225331352663019206},
    }};
    for (double alpha : {1 - 1e-8, 1 + 1e-8}) {
        const StableRand<double> farTail(alpha, 0.5, 1, 0, Parameterization::S0);
        for (const auto &[x, density, cdf] : farReferences) {
            CHECK_THAT(farTail.f(x), WithinAbs(density, 1e-12));
            CHECK_THAT(farTail.F(x), WithinAbs(cdf, 1e-11));
        }
    }
}

TEST_CASE("alpha-one stable tail expansions preserve scale and skew direction", "[stable]")
{
    // Generated by tools/refdata/distributions.py with mpmath at 40 decimal digits.
    struct Reference
    {
        StableRand<double> distribution;
        double x;
        double density;
        double cdf;
    };
    const std::array references{
        Reference{StableRand<double>(1, 0.5, 2, 1), 10000, 9.561292617879164792984036e-9, 0.9999044440279959911893705},
        Reference{StableRand<double>(1, -0.5, 2, 1), -10000, 9.557467091614584637974787e-9,
                  0.00009553685324468203888854466},
        Reference{StableRand<double>(1, -0.75, 0.5, -2), -1000, 2.804403509017000147108865e-7,
                  0.0002795126110786209635510265},
        Reference{StableRand<double>(1, 0.75, 0.5, -2), 1000, 2.782031258189853700595344e-7,
                  0.9997216046669159827840431},
    };
    for (const auto &[distribution, x, density, cdf] : references) {
        CHECK_THAT(distribution.f(x), WithinAbs(density, 1e-13));
        CHECK_THAT(distribution.F(x), WithinAbs(cdf, 1e-11));
        CHECK_THAT(distribution.S(x), WithinAbs(1 - cdf, 1e-11));
    }

    const StableRand<double> landau(1, 1, 1, 0);
    CHECK(landau.f(-100) >= 0);
    CHECK(landau.F(-100) >= 0);
    CHECK(landau.S(-100) <= 1);
}
