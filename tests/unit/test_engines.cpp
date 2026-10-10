// Seeding and output of the bundled random engines and of the standard uniform variate.

#include <catch2/catch_test_macros.hpp>

#include "RandLib.h"

#include <cmath>
#include <limits>
#include <set>
#include <thread>
#include <vector>

namespace
{

template <class Engine>
std::vector<unsigned long long> stream(unsigned long seed, size_t n = 8)
{
    Engine engine;
    engine.Reseed(seed);
    std::vector<unsigned long long> out(n);
    for (auto &x : out)
        x = engine.Next();
    return out;
}

template <class Engine>
double correlationOfAdjacentSeeds(unsigned long seed)
{
    Engine a, b;
    a.Reseed(seed);
    b.Reseed(seed + 1);
    const int n = 20000;
    double scale = static_cast<double>(a.MaxValue());
    double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
    for (int i = 0; i < n; ++i) {
        double x = a.Next() / scale, y = b.Next() / scale;
        sa += x;
        sb += y;
        saa += x * x;
        sbb += y * y;
        sab += x * y;
    }
    double cov = sab / n - sa / n * sb / n;
    return cov / std::sqrt((saa / n - sa / n * sa / n) * (sbb / n - sb / n * sb / n));
}

template <class Engine>
void checkSeeding()
{
    CHECK(stream<Engine>(7) == stream<Engine>(7));
    CHECK(stream<Engine>(7) != stream<Engine>(8));

    std::set<std::vector<unsigned long long>> seen;
    for (unsigned long seed = 0; seed < 2000; ++seed)
        seen.insert(stream<Engine>(seed, 2));
    CHECK(seen.size() == 2000);

    if (sizeof(unsigned long) >= 8) {
        unsigned long high = static_cast<unsigned long>(1ULL << 32);
        CHECK(stream<Engine>(7) != stream<Engine>(7 + high));
    }

    for (unsigned long seed : {0UL, 1UL, 12345UL})
        CHECK(std::fabs(correlationOfAdjacentSeeds<Engine>(seed)) < 0.05);

    Engine first, second;
    CHECK(first.Next() != second.Next());
}

}

TEST_CASE("JKiss engine seeding", "[engines]")
{
    checkSeeding<JKissRandEngine>();
}

TEST_CASE("JLKiss64 engine seeding", "[engines]")
{
    checkSeeding<JLKiss64RandEngine>();
}

TEST_CASE("PCG engine seeding", "[engines]")
{
    checkSeeding<PCGRandEngine>();
}

TEST_CASE("PCG engine follows the reference seeding", "[engines]")
{
    /// pcg32_srandom_r(&rng, 42u, 54u) followed by pcg32_random_r, from pcg32-demo.c
    PCGRandEngine engine;
    engine.Reseed(42, 54);
    const unsigned long long expected[] = {0xa15c02b7, 0x7b47f409, 0xba1d3330, 0x83d2f293, 0xbfa4784b, 0xcbed606e};
    for (unsigned long long value : expected)
        CHECK(engine.Next() == value);
}

TEST_CASE("Default-seeded engines differ across threads", "[engines]")
{
    unsigned long long a = 0, b = 0;
    std::thread first([&a] { JKissRandEngine e; a = e.Next(); });
    std::thread second([&b] { JKissRandEngine e; b = e.Next(); });
    first.join();
    second.join();
    CHECK(a != b);
}

TEST_CASE("Standard uniform variates lie in the open unit interval", "[engines]")
{
    RandGenerator generator;
    generator.Reseed(1);
    for (int i = 0; i < 100000; ++i) {
        double u = UniformRand<double>::StandardVariate(generator);
        REQUIRE(u > 0.0);
        REQUIRE(u < 1.0);
        float uf = UniformRand<float>::StandardVariate(generator);
        REQUIRE(uf > 0.0f);
        REQUIRE(uf < 1.0f);
        long double ul = UniformRand<long double>::StandardVariate(generator);
        REQUIRE(ul > 0.0L);
        REQUIRE(ul < 1.0L);
    }
}

TEST_CASE("Standard uniform variates have full double resolution", "[engines]")
{
    RandGenerator generator;
    generator.Reseed(2);
    int onCoarseGrid = 0;
    for (int i = 0; i < 1000; ++i) {
        double scaled = UniformRand<double>::StandardVariate(generator) * 4294967296.0;
        if (scaled - std::floor(scaled) == 0.5)
            ++onCoarseGrid;
    }
    CHECK(onCoarseGrid < 10);
}

TEST_CASE("Closed and half-open uniform variates stay in range", "[engines]")
{
    RandGenerator generator;
    generator.Reseed(3);
    for (int i = 0; i < 100000; ++i) {
        double closed = UniformRand<double>::StandardVariateClosed(generator);
        REQUIRE(closed >= 0.0);
        REQUIRE(closed <= 1.0);
        double halfOpen = UniformRand<double>::StandardVariateHalfClosed(generator);
        REQUIRE(halfOpen >= 0.0);
        REQUIRE(halfOpen < 1.0);
    }
}
