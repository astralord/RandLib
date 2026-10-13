// Quadrature, root finding and optimization in randlib::math.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <cmath>
#include <numbers>
#include <span>
#include <vector>

using namespace randlib;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Adaptive quadrature on finite and infinite intervals", "[numerics]")
{
    CHECK_THAT(math::integrate([] (double x) { return std::exp(x); }, 0, 1), WithinRel(std::numbers::e - 1, 1e-13));
    CHECK_THAT(math::integrate([] (double x) { return 1 / std::sqrt(x); }, 0, 1, 1e-12, 4096), WithinRel(2.0, 1e-10));
    CHECK_THAT(math::integrate([] (double x) { return std::exp(-x * x); }, -INFINITY, INFINITY), WithinRel(std::sqrt(std::numbers::pi), 1e-12));
    CHECK_THAT(math::integrate([] (double x) { return std::exp(-x); }, 2, INFINITY), WithinRel(std::exp(-2.0), 1e-12));
    CHECK_THAT(math::integrate([] (double x) { return 1 / (1 + x * x); }, -INFINITY, 0), WithinRel(std::numbers::pi / 2, 1e-12));
    CHECK(math::integrate([] (double x) { return x; }, 1, 0) == -math::integrate([] (double x) { return x; }, 0, 1));
    CHECK_THAT(static_cast<double>(math::integrate<long double>([] (long double x) { return x * x; }, 0, 3)), WithinRel(9.0, 1e-15));
}

TEST_CASE("TOMS 748 finds bracketed roots to full precision", "[numerics]")
{
    double root = math::findRoot([] (double x) { return std::cos(x) - x; }, 0, 1);
    CHECK_THAT(root, WithinRel(0.73908513321516064166, 4e-16));
    double cubic = math::findRoot([] (double x) { return x * x * x - 2; }, 0, 10);
    CHECK_THAT(cubic, WithinRel(std::cbrt(2.0), 4e-16));
    /// a steep function with a root near one end
    double steep = math::findRoot([] (double x) { return std::exp(20 * x) - 1e-3; }, -5, 5);
    CHECK_THAT(steep, WithinRel(std::log(1e-3) / 20, 1e-14));
    /// a discontinuous function: the bracket closes on the jump
    double jump = math::findRoot([] (double x) { return x < 0.3 ? -1.0 : 1.0; }, 0, 1);
    CHECK_THAT(jump, WithinAbs(0.3, 1e-15));
    CHECK_THROWS_AS(math::findRoot([] (double x) { return x * x + 1; }, -1, 1), std::invalid_argument);
}

TEST_CASE("Multivariate minimization", "[numerics]")
{
    auto rosenbrock = [] (std::span<const double> x) {
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

    auto quadratic = [] (std::span<const double> x) {
        double sum = 0;
        for (size_t i = 0; i < x.size(); ++i)
            sum += (i + 1) * (x[i] - i) * (x[i] - i);
        return sum;
    };
    auto result = math::bfgs(quadratic, std::vector<double>(5, 10.0));
    for (size_t i = 0; i < 5; ++i)
        CHECK_THAT(result.x[i], WithinAbs(static_cast<double>(i), 1e-6));
}
