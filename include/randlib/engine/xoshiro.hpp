#pragma once

#include <randlib/engine/splitmix64.hpp>

#include <array>
#include <bit>
#include <concepts>
#include <cstdint>

namespace randlib
{

/**
 * @brief The Xoshiro256pp class <BR>
 * xoshiro256++ 1.0 by David Blackman and Sebastiano Vigna: 256 bits of state, period 2^256 - 1.
 * jump() advances the state by 2^128 steps and long_jump() by 2^192 steps, which gives
 * non-overlapping streams for parallel simulations (see make_streams).
 */
class Xoshiro256pp
{
public:
    using result_type = std::uint64_t;
    static constexpr std::uint64_t default_seed = 0x853C49E6748FEA9BULL;

    constexpr explicit Xoshiro256pp(std::uint64_t value = default_seed) noexcept { seed(value); }

    template <class SeedSeq>
        requires(!std::convertible_to<SeedSeq &, std::uint64_t>)
    explicit Xoshiro256pp(SeedSeq &sequence)
    {
        seed(sequence);
    }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    /// Expands the seed with SplitMix64, as recommended by the authors
    constexpr void seed(std::uint64_t value) noexcept
    {
        SplitMix64 expand(value);
        for (std::uint64_t &word : s)
            word = expand();
    }

    template <class SeedSeq>
        requires(!std::convertible_to<SeedSeq &, std::uint64_t>)
    void seed(SeedSeq &sequence)
    {
        std::array<std::uint32_t, 8> words{};
        sequence.generate(words.begin(), words.end());
        for (std::size_t i = 0; i < s.size(); ++i)
            s[i] = (static_cast<std::uint64_t>(words[2 * i]) << 32) | words[2 * i + 1];
        if (s[0] == 0 && s[1] == 0 && s[2] == 0 && s[3] == 0)
            s[0] = 1;
    }

    constexpr result_type operator()() noexcept
    {
        const std::uint64_t result = std::rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = std::rotl(s[3], 45);
        return result;
    }

    constexpr void discard(unsigned long long n) noexcept
    {
        while (n-- > 0)
            (*this)();
    }

    /// Advances the state by 2^128 steps
    constexpr void jump() noexcept
    {
        polynomialJump({0x180EC6D33CFD0ABAULL, 0xD5A61266F0C9392CULL, 0xA9582618E03FC9AAULL, 0x39ABDC4529B1661CULL});
    }

    /// Advances the state by 2^192 steps
    constexpr void long_jump() noexcept
    {
        polynomialJump({0x76E15D3EFEFDCBBFULL, 0xC5004E441C522FB3ULL, 0x77710069854EE241ULL, 0x39109BB02ACBE635ULL});
    }

    friend constexpr bool operator==(const Xoshiro256pp &, const Xoshiro256pp &) = default;

private:
    std::array<std::uint64_t, 4> s{};

    constexpr void polynomialJump(const std::array<std::uint64_t, 4> &polynomial) noexcept
    {
        std::array<std::uint64_t, 4> next{};
        for (std::uint64_t word : polynomial) {
            for (int b = 0; b < 64; ++b) {
                if (word & (std::uint64_t{1} << b)) {
                    for (std::size_t i = 0; i < s.size(); ++i)
                        next[i] ^= s[i];
                }
                (*this)();
            }
        }
        s = next;
    }
};

} // namespace randlib
