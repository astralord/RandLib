#pragma once

#include <randlib/engine/pcg.hpp>
#include <randlib/engine/splitmix64.hpp>
#include <randlib/engine/xoshiro.hpp>

#include <chrono>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <thread>
#include <vector>

namespace randlib
{

/// The engine used by default, including by the per-instance engines of the distributions
using DefaultEngine = Xoshiro256pp;

/**
 * @fn random_seed
 * @return non-deterministic 64-bit seed: std::random_device mixed with a high-resolution clock and
 * the thread id through SplitMix64
 */
inline std::uint64_t random_seed()
{
    std::uint64_t entropy = 0;
    try {
        std::random_device device;
        entropy = (static_cast<std::uint64_t>(device()) << 32) ^ device();
    } catch (const std::exception &) {
    }
    entropy ^= static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    entropy ^= static_cast<std::uint64_t>(std::hash<std::thread::id>()(std::this_thread::get_id())) << 1;
    return SplitMix64(entropy)();
}

/**
 * @fn seed_from
 * @param sequence seed sequence that fills the whole state of the engine
 * @return engine seeded from the sequence
 */
template <class Engine = DefaultEngine>
Engine seed_from(std::seed_seq &sequence)
{
    return Engine(sequence);
}

/**
 * @fn make_streams
 * Engines for parallel simulations, one per thread: consecutive jump() calls for engines that
 * have them (non-overlapping by construction), distinct PCG streams for Pcg32 and Pcg64, and
 * seeds derived with SplitMix64 for any other engine.
 * @param seed common seed
 * @param count number of engines
 */
template <class Engine = DefaultEngine>
std::vector<Engine> make_streams(std::uint64_t seed, std::size_t count)
{
    std::vector<Engine> streams;
    streams.reserve(count);
    if constexpr (requires(Engine e) { e.jump(); }) {
        Engine engine(seed);
        for (std::size_t i = 0; i < count; ++i) {
            streams.push_back(engine);
            engine.jump();
        }
    } else if constexpr (std::constructible_from<Engine, std::uint64_t, std::uint64_t>) {
        for (std::size_t i = 0; i < count; ++i)
            streams.emplace_back(seed, static_cast<std::uint64_t>(i));
    } else {
        SplitMix64 seeds(seed);
        for (std::size_t i = 0; i < count; ++i)
            streams.emplace_back(static_cast<typename Engine::result_type>(seeds()));
    }
    return streams;
}

} // namespace randlib
