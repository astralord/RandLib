// The same seed gives the same sequence, independently of other instances and of call
// history: drawing from another instance, or asking for quantiles, must not change it.

#include <catch2/catch_test_macros.hpp>

#include "registry.h"

#include <vector>

using namespace randlib_test;

namespace
{

constexpr unsigned long SEED = 42;
constexpr int DRAWS = 300;

template <class Make>
void checkReproducible(const Case &c, Make make)
{
    INFO(c.id);
    decltype(make()) a, b, other;
    try {
        a = make();
        b = make();
        other = make();
    } catch (const std::exception &e) {
        FAIL_CHECK(c.id + ": construction failed: " + e.what());
        return;
    }

    a->Reseed(SEED);
    std::vector<double> first;
    for (int i = 0; i < DRAWS; ++i)
        first.push_back(a->Variate());

    other->Reseed(SEED + 1);
    for (int i = 0; i < DRAWS; ++i)
        other->Variate();

    b->Reseed(SEED);
    std::vector<double> second;
    for (int i = 0; i < DRAWS; ++i)
        second.push_back(b->Variate());
    CHECK(first == second);

    b->Reseed(SEED);
    try {
        b->Median();
        b->Quantile(0.3);
        b->Quantile1m(0.2);
    } catch (const randlib::not_implemented &) {
    }
    std::vector<double> third;
    for (int i = 0; i < DRAWS; ++i)
        third.push_back(b->Variate());
    CHECK(first == third);
}

}

TEST_CASE("Same seed reproduces the same continuous sequence", "[reproducibility]")
{
    for (const Case *c : casesOfKind(Kind::Continuous))
        checkReproducible(*c, [c] { return c->makeContinuous(); });
}

TEST_CASE("Same seed reproduces the same discrete sequence", "[reproducibility]")
{
    for (const Case *c : casesOfKind(Kind::Discrete))
        checkReproducible(*c, [c] { return c->makeDiscrete(); });
}

TEST_CASE("Same seed reproduces the same singular sequence", "[reproducibility]")
{
    for (const Case *c : casesOfKind(Kind::Singular))
        checkReproducible(*c, [c] { return c->makeSingular(); });
}

TEST_CASE("Bulk sampling is reproducible", "[reproducibility]")
{
    for (const Case *c : casesOfKind(Kind::Continuous)) {
        INFO(c->id);
        std::unique_ptr<ContinuousDistribution<double>> a, other;
        try {
            a = c->makeContinuous();
            other = c->makeContinuous();
        } catch (const std::exception &e) {
            FAIL_CHECK(c->id + ": construction failed: " + e.what());
            continue;
        }
        std::vector<double> first(DRAWS), second(DRAWS);
        a->Reseed(SEED);
        a->Sample(first);
        other->Reseed(SEED);
        other->Sample(second);
        CHECK(first == second);
    }
    for (const Case *c : casesOfKind(Kind::Discrete)) {
        INFO(c->id);
        auto a = c->makeDiscrete(), other = c->makeDiscrete();
        std::vector<int> first(DRAWS), second(DRAWS);
        a->Reseed(SEED);
        a->Sample(first);
        other->Reseed(SEED);
        other->Sample(second);
        CHECK(first == second);
    }
}
