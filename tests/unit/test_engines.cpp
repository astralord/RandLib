// Bundled random engines, seeding helpers and the conversion of engine output into uniform variates.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

using namespace randlib;
using Catch::Matchers::WithinULP;

static_assert(std::uniform_random_bit_generator<SplitMix64>);
static_assert(std::uniform_random_bit_generator<Pcg32>);
static_assert(std::uniform_random_bit_generator<Pcg64>);
static_assert(std::uniform_random_bit_generator<Xoshiro256pp>);
static_assert(std::uniform_random_bit_generator<legacy::JKissEngine>);
static_assert(std::uniform_random_bit_generator<legacy::JLKiss64Engine>);
static_assert(std::uniform_random_bit_generator<detail::EngineRef>);

namespace
{

template <class Engine>
std::vector<std::uint64_t> stream(std::uint64_t seed, size_t n = 8)
{
    Engine engine(seed);
    std::vector<std::uint64_t> out(n);
    for (auto &x : out)
        x = engine();
    return out;
}

template <class Engine>
void checkSeeding()
{
    CHECK(stream<Engine>(7) == stream<Engine>(7));
    CHECK(stream<Engine>(7) != stream<Engine>(8));
    CHECK(stream<Engine>(7) != stream<Engine>(7 + (std::uint64_t{1} << 32)));

    std::set<std::vector<std::uint64_t>> seen;
    for (std::uint64_t seed = 0; seed < 2000; ++seed)
        seen.insert(stream<Engine>(seed, 2));
    CHECK(seen.size() == 2000);

    Engine reseeded(1);
    reseeded.seed(7);
    CHECK(reseeded == Engine(7));
}

} // namespace

TEST_CASE("Engines are seeded from all bits of the seed", "[engines]")
{
    checkSeeding<Xoshiro256pp>();
    checkSeeding<Pcg32>();
    checkSeeding<Pcg64>();
    checkSeeding<legacy::JKissEngine>();
    checkSeeding<legacy::JLKiss64Engine>();
}

TEST_CASE("PCG engines follow the reference implementations", "[engines]")
{
    /// pcg32_srandom_r(&rng, 42u, 54u) followed by pcg32_random_r, from pcg32-demo.c
    Pcg32 pcg32(42, 54);
    for (std::uint32_t value : {0xa15c02b7u, 0x7b47f409u, 0xba1d3330u, 0x83d2f293u, 0xbfa4784bu, 0xcbed606eu})
        CHECK(pcg32() == value);

    /// pcg64 seeded with (42, 54), checked against NumPy's PCG64 set to the same state and increment
    Pcg64 pcg64(42, 54);
    for (std::uint64_t value : {9705778491962043240ULL, 1370407407632858425ULL, 11774395822783136600ULL})
        CHECK(pcg64() == value);
}

TEST_CASE("Xoshiro256++ follows the reference implementation", "[engines]")
{
    /// state expanded from seed 1 by SplitMix64, as the authors advise
    Xoshiro256pp engine(1);
    CHECK(engine() == 14971601782005023387ULL);

    /// jumps are polynomials in the transition, so they commute with it
    Xoshiro256pp jumpThenStep(5), stepThenJump(5);
    jumpThenStep.jump();
    jumpThenStep();
    stepThenJump();
    stepThenJump.jump();
    CHECK(jumpThenStep == stepThenJump);

    Xoshiro256pp jumped(5), longJumped(5);
    jumped.jump();
    longJumped.long_jump();
    CHECK(jumped != Xoshiro256pp(5));
    CHECK(longJumped != jumped);
}

TEST_CASE("Legacy engines reproduce the raw output of RandLib 1.2", "[engines]")
{
    legacy::JKissEngine jkiss(7);
    for (std::uint32_t value : {0x41e0ca2au, 0xe36ed7bu, 0xa3e4ae72u, 0xe5fb9220u})
        CHECK(jkiss() == value);
    legacy::JLKiss64Engine jlkiss(7);
    for (std::uint64_t value :
         {0x3522a7f986f0823fULL, 0x1d3111d61d5c3260ULL, 0x5ef9e390aef7967eULL, 0x208fc63648fab68cULL})
        CHECK(jlkiss() == value);
}

TEST_CASE("Seeding helpers", "[engines]")
{
    CHECK(random_seed() != random_seed());

    std::seed_seq sequence{1, 2, 3};
    auto first = seed_from(sequence);
    auto second = seed_from(sequence);
    CHECK(first == second);
    CHECK(first != DefaultEngine());

    auto streams = make_streams(11, 3);
    Xoshiro256pp expected(11);
    for (auto &engine : streams) {
        CHECK(engine == expected);
        expected.jump();
    }

    auto pcgStreams = make_streams<Pcg64>(11, 2);
    CHECK(pcgStreams[0] == Pcg64(11, 0));
    CHECK(pcgStreams[1] == Pcg64(11, 1));

    auto mtStreams = make_streams<std::mt19937_64>(11, 2);
    CHECK(mtStreams[0]() != mtStreams[1]());
}

TEST_CASE("Random bits are uniform for engines of any range", "[engines]")
{
    /// minstd_rand covers [1, 2^31 - 2], which is not a power of two
    std::minstd_rand engine(3);
    std::array<int, 52> ones{};
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        std::uint64_t bits = detail::randomBits(engine, 52);
        REQUIRE(bits >> 52 == 0);
        for (int b = 0; b < 52; ++b)
            ones[b] += (bits >> b) & 1;
    }
    for (int count : ones)
        CHECK(std::fabs(count - 0.5 * n) < 5 * std::sqrt(0.25 * n));
}

TEST_CASE("Uniform integers have no modulo bias", "[engines]")
{
    Pcg32 engine(4);
    const std::uint64_t n = 6;
    std::array<int, n> counts{};
    for (int i = 0; i < 60000; ++i)
        ++counts[detail::uniformBelow(engine, n)];
    double chiSquared = 0;
    for (int count : counts)
        chiSquared += (count - 10000.0) * (count - 10000.0) / 10000.0;
    CHECK(chiSquared < 25); /// p < 1e-4 for 5 degrees of freedom
}

TEST_CASE("Standard uniform variates lie in the open unit interval", "[engines]")
{
    std::mt19937 mt(1);
    Pcg64 pcg(1);
    for (int i = 0; i < 100000; ++i) {
        double u = detail::uniform01<double>(mt);
        REQUIRE(u > 0.0);
        REQUIRE(u < 1.0);
        float uf = detail::uniform01<float>(pcg);
        REQUIRE(uf > 0.0f);
        REQUIRE(uf < 1.0f);
        long double ul = detail::uniform01<long double>(mt);
        REQUIRE(ul > 0.0L);
        REQUIRE(ul < 1.0L);
    }
}

TEST_CASE("Standard uniform variates have full double resolution", "[engines]")
{
    /// 32-bit engines need two outputs per double
    std::mt19937 engine(2);
    int onCoarseGrid = 0;
    for (int i = 0; i < 1000; ++i) {
        double scaled = detail::uniform01<double>(engine) * 4294967296.0;
        if (scaled - std::floor(scaled) == 0.5)
            ++onCoarseGrid;
    }
    CHECK(onCoarseGrid < 10);
}

TEST_CASE("Closed and half-open uniform variates stay in range", "[engines]")
{
    Xoshiro256pp engine(3);
    for (int i = 0; i < 100000; ++i) {
        double closed = detail::uniform01Closed<double>(engine);
        REQUIRE(closed >= 0.0);
        REQUIRE(closed <= 1.0);
        double halfOpen = detail::uniform01HalfClosed<double>(engine);
        REQUIRE(halfOpen >= 0.0);
        REQUIRE(halfOpen < 1.0);
    }
}

TEST_CASE("Distributions sample with any engine", "[engines]")
{
    NormalRand<double> X(1, 4);
    std::mt19937 mt(5);
    std::minstd_rand minstd(5);
    double sum = 0;
    for (int i = 0; i < 1000; ++i)
        sum += X(mt) + X(minstd);
    CHECK(std::fabs(sum / 2000 - 1.0) < 0.3);

    /// with a 64-bit engine, the virtual interface consumes the engine like direct calls
    Pcg64 direct(6), bulk(6), throughBase(6);
    std::vector<double> sample(16);
    X.Sample(bulk, sample);
    const ContinuousDistribution<double> &base = X;
    for (double value : sample) {
        CHECK_THAT(X(direct), WithinULP(value, 4));
        CHECK_THAT(base(throughBase), WithinULP(value, 4));
    }
    auto directNext = direct();
    auto bulkNext = bulk();
    auto throughBaseNext = throughBase();
    CHECK(directNext == bulkNext);
    CHECK(throughBaseNext == bulkNext);
}

TEST_CASE("Default-seeded distributions differ", "[engines]")
{
    NormalRand<double> first, second;
    CHECK(first.Variate() != second.Variate());
}
