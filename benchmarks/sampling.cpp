#include <benchmark/benchmark.h>

#include <randlib/randlib.hpp>

#include <cstdint>
#include <random>

namespace
{

template <typename Distribution>
void sample(benchmark::State &state, Distribution distribution)
{
    std::mt19937_64 engine(0x5eed);
    for (auto _ : state)
        benchmark::DoNotOptimize(distribution(engine));
    state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
}

void randlibNormal(benchmark::State &state)
{
    sample(state, randlib::NormalRand<double>(1.5, 4.0));
}

void standardNormal(benchmark::State &state)
{
    sample(state, std::normal_distribution<double>(1.5, 2.0));
}

void randlibExponential(benchmark::State &state)
{
    sample(state, randlib::ExponentialRand<double>(0.75));
}

void standardExponential(benchmark::State &state)
{
    sample(state, std::exponential_distribution<double>(0.75));
}

void randlibGamma(benchmark::State &state)
{
    sample(state, randlib::GammaRand<double>(2.5, 0.75));
}

void standardGamma(benchmark::State &state)
{
    sample(state, std::gamma_distribution<double>(2.5, 1.0 / 0.75));
}

void randlibPoisson(benchmark::State &state)
{
    sample(state, randlib::PoissonRand<int>(12.0));
}

void standardPoisson(benchmark::State &state)
{
    sample(state, std::poisson_distribution<int>(12.0));
}

void randlibStable(benchmark::State &state)
{
    sample(state, randlib::StableRand<double>(1.7, 0.5, 2.0, 1.0));
}

BENCHMARK(randlibNormal);
BENCHMARK(standardNormal);
BENCHMARK(randlibExponential);
BENCHMARK(standardExponential);
BENCHMARK(randlibGamma);
BENCHMARK(standardGamma);
BENCHMARK(randlibPoisson);
BENCHMARK(standardPoisson);
BENCHMARK(randlibStable);

} // namespace
