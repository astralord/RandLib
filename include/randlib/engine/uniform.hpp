#pragma once

#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>

namespace randlib::detail
{

/// Number of uniformly distributed bits in one output of G: floor(log2(max - min + 1))
template <std::uniform_random_bit_generator G>
constexpr int engineBits() noexcept
{
    static_assert(sizeof(typename G::result_type) <= sizeof(std::uint64_t),
                  "engines wider than 64 bits are not supported");
    constexpr std::uint64_t range = static_cast<std::uint64_t>(G::max() - G::min());
    if constexpr (range == std::numeric_limits<std::uint64_t>::max())
        return 64;
    else
        return std::bit_width(range + 1) - 1;
}

/**
 * @param g any uniform random bit generator, including those whose range is not a power of two
 * @param count number of bits, from 1 to 64
 * @return count independent uniform random bits, taken from the high end of the engine outputs
 */
template <std::uniform_random_bit_generator G>
std::uint64_t randomBits(G &g, int count)
{
    constexpr int width = engineBits<G>();
    constexpr std::uint64_t range = static_cast<std::uint64_t>(G::max() - G::min());
    constexpr bool exact = width == 64 || range + 1 == (std::uint64_t{1} << (width % 64));
    if constexpr (width == 64) {
        std::uint64_t v = static_cast<std::uint64_t>(g() - G::min());
        return count == 64 ? v : v >> (64 - count);
    } else {
        std::uint64_t bits = 0;
        int filled = 0;
        while (filled < count) {
            std::uint64_t v = static_cast<std::uint64_t>(g() - G::min());
            if constexpr (!exact) {
                if (v >> width)
                    continue;
            }
            bits = (bits << width) | v;
            filled += width;
        }
        int available = filled < 64 ? filled : 64;
        return count == 64 ? bits : bits >> (available - count);
    }
}

/**
 * @return (k + 1/2) / 2^b for b = digits - 1 random bits k: uniform on (0, 1) with the full
 * resolution of Real, never 0 and never 1
 */
template <std::floating_point Real, std::uniform_random_bit_generator G>
Real uniform01(G &g)
{
    constexpr int b = std::numeric_limits<Real>::digits - 1;
    if constexpr (b <= 64) {
        Real k = static_cast<Real>(randomBits(g, b));
        return std::ldexp(k + static_cast<Real>(0.5), -b);
    } else {
        Real high = static_cast<Real>(randomBits(g, 64));
        Real low = static_cast<Real>(randomBits(g, b - 64));
        return std::ldexp(std::ldexp(high, b - 64) + low + static_cast<Real>(0.5), -b);
    }
}

/// k / (2^b - 1) for b = min(digits, 64) random bits: uniform on [0, 1]
template <std::floating_point Real, std::uniform_random_bit_generator G>
Real uniform01Closed(G &g)
{
    constexpr int b = std::numeric_limits<Real>::digits < 64 ? std::numeric_limits<Real>::digits : 64;
    Real k = static_cast<Real>(randomBits(g, b));
    return k / (std::ldexp(static_cast<Real>(1), b) - 1);
}

/// k / 2^b for b = min(digits, 64) random bits: uniform on [0, 1)
template <std::floating_point Real, std::uniform_random_bit_generator G>
Real uniform01HalfClosed(G &g)
{
    constexpr int b = std::numeric_limits<Real>::digits < 64 ? std::numeric_limits<Real>::digits : 64;
    Real k = static_cast<Real>(randomBits(g, b));
    return std::ldexp(k, -b);
}

/**
 * @param g uniform random bit generator
 * @param n positive number of outcomes
 * @return uniform integer in [0, n), without modulo bias: 64-bit values below 2^64 mod n are rejected
 */
template <std::uniform_random_bit_generator G>
std::uint64_t uniformBelow(G &g, std::uint64_t n)
{
    const std::uint64_t threshold = (0 - n) % n;
    std::uint64_t v = randomBits(g, 64);
    while (v < threshold)
        v = randomBits(g, 64);
    return v % n;
}

/**
 * @brief The EngineRef class <BR>
 * Non-owning, type-erased reference to any uniform random bit generator. Each call returns 64 random
 * bits of the referenced engine, so the virtual sampling interface of the distribution classes can
 * accept engines of every type. Engines with 64-bit output are consumed exactly as by direct calls.
 */
class EngineRef
{
public:
    using result_type = std::uint64_t;

    template <std::uniform_random_bit_generator G>
        requires(!std::same_as<std::remove_cv_t<G>, EngineRef>)
    EngineRef(G &engine) noexcept : target(&engine), next([](void *e) { return randomBits(*static_cast<G *>(e), 64); })
    {
    }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    result_type operator()() { return next(target); }

private:
    void *target;
    result_type (*next)(void *);
};

} // namespace randlib::detail
