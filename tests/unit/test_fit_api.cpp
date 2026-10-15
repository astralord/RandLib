#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <cmath>
#include <vector>

using namespace randlib;
using Catch::Matchers::WithinAbs;

TEST_CASE("Unified fitting reports likelihood and information criteria", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3, 4};
    auto result = fit<NormalRand<double>>(sample, Method::MLE);

    CHECK(result.converged);
    CHECK(result.iterations == 0);
    CHECK(result.method == Method::MLE);
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(2.5, 1e-14));
    CHECK_THAT(static_cast<double>(result.distribution.Variance()), WithinAbs(1.25, 1e-14));
    CHECK(result.parameter_names == std::vector<std::string>{"location", "variance"});
    REQUIRE(result.parameters.size() == 2);
    CHECK_THAT(result.parameters[0], WithinAbs(2.5, 1e-14));
    CHECK_THAT(result.parameters[1], WithinAbs(1.25, 1e-14));
    CHECK(result.std_errors.size() == 2);
    CHECK(std::isfinite(result.std_errors[0]));
    CHECK(std::isfinite(result.std_errors[1]));
    CHECK_THAT(result.aic, WithinAbs(4 - 2 * result.log_likelihood, 1e-14));
    CHECK_THAT(result.bic, WithinAbs(2 * std::log(4.0) - 2 * result.log_likelihood, 1e-14));
}

TEST_CASE("Unified fitting supports fixed parameters", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3, 4};
    FitOptions options;
    options.fixed = {{"variance", 2.0}};
    auto result = fit<NormalRand<double>>(sample, options);

    CHECK(result.converged);
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(2.5, 1e-7));
    CHECK_THAT(static_cast<double>(result.distribution.Variance()), WithinAbs(2.0, 1e-14));
    CHECK(result.fixed == std::vector<bool>{false, true});
    CHECK(result.std_errors[1] == 0);
    CHECK_THAT(result.aic, WithinAbs(2 - 2 * result.log_likelihood, 1e-14));
}

TEST_CASE("Unified fitting supports discrete distributions", "[fit]")
{
    const std::vector<int> sample = {0, 1, 2, 3};
    auto result = fit<PoissonRand<int>>(sample);

    CHECK(result.converged);
    CHECK_THAT(result.distribution.GetRate(), WithinAbs(1.5, 1e-14));
    CHECK(result.parameter_names == std::vector<std::string>{"rate"});
    CHECK(result.std_errors[0] > 0);
}

TEST_CASE("Unified fitting rejects invalid options", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3};
    FitOptions options;
    options.fixed = {{"unknown", 1.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), std::invalid_argument);
    options.fixed = {{"variance", -1.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), std::invalid_argument);
    options.fixed = {{"variance", 0.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), std::invalid_argument);
    options.method = Method::UMVU;
    options.fixed = {{"variance", 1.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), fit_error);
    CHECK_THROWS_AS(fit<NormalRand<double>>(std::span<const double>{}), fit_error);
}
