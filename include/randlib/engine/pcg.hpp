#pragma once

#include <bit>
#include <concepts>
#include <cstdint>

namespace randlib
{

/**
 * @brief The Pcg32 class <BR>
 * PCG-XSH-RR with 64 bits of state and 32-bit output (Melissa O'Neill, https://www.pcg-random.org),
 * seeded exactly as the reference pcg32_srandom_r(seed, stream).
 */
class Pcg32
{
public:
    using result_type = std::uint32_t;
    static constexpr std::uint64_t default_seed = 0x853C49E6748FEA9BULL;
    static constexpr std::uint64_t default_stream = 0xDA3E39CB94B95BDBULL;

    constexpr explicit Pcg32(std::uint64_t value = default_seed, std::uint64_t stream = default_stream) noexcept
    {
        seed(value, stream);
    }

    template <class SeedSeq>
        requires (!std::convertible_to<SeedSeq &, std::uint64_t>)
    explicit Pcg32(SeedSeq &sequence) { seed(sequence); }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    constexpr void seed(std::uint64_t value, std::uint64_t stream = default_stream) noexcept
    {
        state = 0;
        increment = (stream << 1) | 1;
        step();
        state += value;
        step();
    }

    template <class SeedSeq>
        requires (!std::convertible_to<SeedSeq &, std::uint64_t>)
    void seed(SeedSeq &sequence)
    {
        std::uint32_t words[4] = {};
        sequence.generate(words, words + 4);
        seed((static_cast<std::uint64_t>(words[0]) << 32) | words[1], (static_cast<std::uint64_t>(words[2]) << 32) | words[3]);
    }

    constexpr result_type operator()() noexcept
    {
        std::uint64_t previous = state;
        step();
        auto xorShifted = static_cast<std::uint32_t>(((previous >> 18) ^ previous) >> 27);
        return std::rotr(xorShifted, static_cast<int>(previous >> 59));
    }

    constexpr void discard(unsigned long long n) noexcept
    {
        while (n-- > 0)
            step();
    }

    friend constexpr bool operator==(const Pcg32 &, const Pcg32 &) = default;

private:
    static constexpr std::uint64_t MULTIPLIER = 6364136223846793005ULL;
    std::uint64_t state = 0;
    std::uint64_t increment = 1;

    constexpr void step() noexcept { state = state * MULTIPLIER + increment; }
};

/**
 * @brief The Pcg64 class <BR>
 * PCG-XSL-RR with 128 bits of state and 64-bit output, the pcg64 engine of the reference C++
 * library and of NumPy, seeded as pcg64(seed, stream).
 */
class Pcg64
{
public:
    using result_type = std::uint64_t;
    static constexpr std::uint64_t default_seed = 0xCAFEF00DD15EA5E5ULL;
    static constexpr std::uint64_t default_stream = 0;

    constexpr explicit Pcg64(std::uint64_t value = default_seed, std::uint64_t stream = default_stream) noexcept
    {
        seed(value, stream);
    }

    template <class SeedSeq>
        requires (!std::convertible_to<SeedSeq &, std::uint64_t>)
    explicit Pcg64(SeedSeq &sequence) { seed(sequence); }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    /// pcg_engines::setseq_xsl_rr_128_64(seed, stream) with 64-bit seed and stream
    constexpr void seed(std::uint64_t value, std::uint64_t stream = default_stream) noexcept
    {
        seed(UInt128{0, value}, UInt128{0, stream});
    }

    template <class SeedSeq>
        requires (!std::convertible_to<SeedSeq &, std::uint64_t>)
    void seed(SeedSeq &sequence)
    {
        std::uint32_t w[8] = {};
        sequence.generate(w, w + 8);
        auto word = [&w](int i) { return (static_cast<std::uint64_t>(w[i]) << 32) | w[i + 1]; };
        seed(UInt128{word(0), word(2)}, UInt128{word(4), word(6)});
    }

    constexpr result_type operator()() noexcept
    {
        step();
        return std::rotr(state.hi ^ state.lo, static_cast<int>(state.hi >> 58));
    }

    constexpr void discard(unsigned long long n) noexcept
    {
        while (n-- > 0)
            step();
    }

    friend constexpr bool operator==(const Pcg64 &, const Pcg64 &) = default;

private:
    struct UInt128
    {
        std::uint64_t hi = 0;
        std::uint64_t lo = 0;

        friend constexpr bool operator==(const UInt128 &, const UInt128 &) = default;

        friend constexpr UInt128 operator+(UInt128 a, UInt128 b) noexcept
        {
            std::uint64_t lowSum = a.lo + b.lo;
            return {a.hi + b.hi + (lowSum < a.lo ? 1 : 0), lowSum};
        }

        friend constexpr UInt128 operator*(UInt128 a, UInt128 b) noexcept
        {
            /// low 128 bits of the product: a.lo * b.lo in full, plus the cross terms shifted by 64
            std::uint64_t a0 = a.lo & 0xFFFFFFFFULL, a1 = a.lo >> 32, b0 = b.lo & 0xFFFFFFFFULL, b1 = b.lo >> 32;
            std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
            std::uint64_t middle = (p00 >> 32) + (p01 & 0xFFFFFFFFULL) + (p10 & 0xFFFFFFFFULL);
            std::uint64_t lowProduct = (middle << 32) | (p00 & 0xFFFFFFFFULL);
            std::uint64_t highProduct = p11 + (p01 >> 32) + (p10 >> 32) + (middle >> 32);
            highProduct += a.hi * b.lo + a.lo * b.hi;
            return {highProduct, lowProduct};
        }
    };

    static constexpr UInt128 MULTIPLIER{0x2360ED051FC65DA4ULL, 0x4385DF649FCCF645ULL};
    UInt128 state{};
    UInt128 increment{0, 1};

    constexpr void seed(UInt128 value, UInt128 stream) noexcept
    {
        state = UInt128{};
        increment = UInt128{(stream.hi << 1) | (stream.lo >> 63), (stream.lo << 1) | 1};
        step();
        state = state + value;
        step();
    }

    constexpr void step() noexcept { state = state * MULTIPLIER + increment; }
};

} // namespace randlib
