#pragma once

#include <cstdint>

namespace randlib
{

/**
 * @brief The SplitMix64 class <BR>
 * Sebastiano Vigna's SplitMix64 generator. It is used to expand a single seed into the
 * state of the other engines and should not be used for simulation on its own.
 */
class SplitMix64
{
public:
    using result_type = std::uint64_t;

    constexpr explicit SplitMix64(std::uint64_t seed = 0) noexcept : state(seed) {}

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    constexpr void seed(std::uint64_t value) noexcept { state = value; }

    constexpr result_type operator()() noexcept
    {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    constexpr void discard(unsigned long long n) noexcept { state += n * 0x9E3779B97F4A7C15ULL; }

    friend constexpr bool operator==(const SplitMix64 &, const SplitMix64 &) = default;

private:
    std::uint64_t state;
};

} // namespace randlib
