#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/continuous/stable.hpp>

#include <array>
#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
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
    CHECK(stable.MinValue() == 1.0);
    CHECK(stable.Name() == "Stable(0.5, 1, 2, 3, S0)");

    stable.SetScale(4.0);
    CHECK(stable.GetLocation() == 3.0);
    CHECK(stable.MinValue() == -1.0);
    stable.SetLocation(7.0);
    CHECK(stable.GetLocation() == 7.0);
    CHECK(stable.MinValue() == 3.0);

    const auto characteristicFunction = stable.CF(0.75);
    stable.SetParameterization(Parameterization::S1);
    CHECK(stable.GetLocation() == 3.0);
    CHECK(stable.CF(0.75) == characteristicFunction);
    CHECK(stable.Name() == "Stable(0.5, 1, 4, 3)");
}
