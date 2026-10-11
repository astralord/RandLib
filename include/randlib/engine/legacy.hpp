#pragma once

#include <randlib/engine/splitmix64.hpp>

#include <cstdint>

/// Engines of RandLib 1.x, kept for comparison with older results. They produce the same raw
/// output as the 1.2 engines for the same seed; sampled values still differ from 1.x, because
/// the samplers changed. New code should use DefaultEngine.
namespace randlib::legacy
{

/**
 * @brief The JKissEngine class <BR>
 * David Jones' JKISS generator with 32-bit output
 */
class JKissEngine
{
public:
    using result_type = std::uint32_t;

    explicit JKissEngine(std::uint64_t value = 0) noexcept { seed(value); }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    void seed(std::uint64_t value) noexcept
    {
        SplitMix64 expand(value);
        X = static_cast<std::uint32_t>(expand());
        do {
            Y = static_cast<std::uint32_t>(expand()); /// xorshift state must not be zero
        } while (Y == 0);
        Z = static_cast<std::uint32_t>(expand());
        C = 1 + static_cast<std::uint32_t>(expand() % 698769067ULL); /// carry in [1, a - 2]
    }

    result_type operator()() noexcept
    {
        std::uint64_t t = 698769069ULL * Z + C;
        X = 69069 * X + 12345;
        Y ^= Y << 13;
        Y ^= Y >> 17;
        Y ^= Y << 5;
        C = static_cast<std::uint32_t>(t >> 32);
        Z = static_cast<std::uint32_t>(t);
        return X + Y + Z;
    }

    friend bool operator==(const JKissEngine &, const JKissEngine &) = default;

private:
    std::uint32_t X = 0, Y = 1, Z = 0, C = 1;
};

/**
 * @brief The JLKiss64Engine class <BR>
 * David Jones' JLKISS64 generator with 64-bit output
 */
class JLKiss64Engine
{
public:
    using result_type = std::uint64_t;

    explicit JLKiss64Engine(std::uint64_t value = 0) noexcept { seed(value); }

    static constexpr result_type min() noexcept { return 0; }
    static constexpr result_type max() noexcept { return ~result_type{0}; }

    void seed(std::uint64_t value) noexcept
    {
        SplitMix64 expand(value);
        X = expand();
        do {
            Y = expand(); /// xorshift state must not be zero
        } while (Y == 0);
        Z1 = static_cast<std::uint32_t>(expand());
        Z2 = static_cast<std::uint32_t>(expand());
        C1 = 1 + static_cast<std::uint32_t>(expand() % 4294584391ULL); /// carries in [1, a - 2]
        C2 = 1 + static_cast<std::uint32_t>(expand() % 4246477507ULL);
    }

    result_type operator()() noexcept
    {
        X = 1490024343005336237ULL * X + 123456789;
        Y ^= Y << 21;
        Y ^= Y >> 17;
        Y ^= Y << 30;
        std::uint64_t t = 4294584393ULL * Z1 + C1;
        C1 = static_cast<std::uint32_t>(t >> 32);
        Z1 = static_cast<std::uint32_t>(t);
        t = 4246477509ULL * Z2 + C2;
        C2 = static_cast<std::uint32_t>(t >> 32);
        Z2 = static_cast<std::uint32_t>(t);
        return X + Y + Z1 + (static_cast<std::uint64_t>(Z2) << 32);
    }

    friend bool operator==(const JLKiss64Engine &, const JLKiss64Engine &) = default;

private:
    std::uint64_t X = 0, Y = 1;
    std::uint32_t Z1 = 0, Z2 = 0, C1 = 1, C2 = 1;
};

} // namespace randlib::legacy
