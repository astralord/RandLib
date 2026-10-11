#pragma once

#include <randlib/core/ziggurat_tables.hpp>
#include <randlib/engine/uniform.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <random>
#include <stdexcept>

namespace randlib::detail
{

/// Upper bound on the iterations of rejection samplers; reaching it means there is a bug
inline constexpr std::size_t MAX_ITER_REJECTION = 1000;

/**
 * @fn standardExponential
 * Ziggurat algorithm of Marsaglia and Tsang with 256 layers
 * @return Exp(1) variate
 */
template <std::floating_point Real, std::uniform_random_bit_generator G>
Real standardExponential(G &g)
{
    const auto &table = exponentialZiggurat;
    for (std::size_t iter = 0; iter <= MAX_ITER_REJECTION; ++iter) {
        int layer = static_cast<int>(randomBits(g, 8));
        Real x = uniform01<Real>(g) * static_cast<Real>(table[layer].second);
        if (x < table[layer + 1].second)
            return x;
        if (layer == 0)
            return static_cast<Real>(table[1].second) + standardExponential<Real>(g);
        Real height = static_cast<Real>(table[layer].first - table[layer - 1].first);
        if (table[layer - 1].first + height * uniform01<Real>(g) < std::exp(-x))
            return x;
    }
    throw std::runtime_error("Exponential distribution: sampling failed");
}

/**
 * @fn standardNormal
 * Ziggurat algorithm of Marsaglia and Tsang with 256 layers and an explicit sign bit
 * @return N(0, 1) variate
 */
template <std::floating_point Real, std::uniform_random_bit_generator G>
Real standardNormal(G &g)
{
    const auto &table = normalZiggurat;
    for (std::size_t iter = 0; iter <= MAX_ITER_REJECTION; ++iter) {
        auto bits = randomBits(g, 9);
        int layer = static_cast<int>(bits >> 1);
        bool positive = bits & 1;
        Real x = uniform01<Real>(g) * static_cast<Real>(table[layer].second);
        if (x < table[layer + 1].second)
            return positive ? x : -x;
        if (layer == 0) {
            /// Marsaglia's method for the tail beyond r: accept x if E > x^2 / 2
            const Real r = static_cast<Real>(table[1].second);
            Real e = 0;
            do {
                x = standardExponential<Real>(g) / r;
                e = standardExponential<Real>(g);
            } while (e + e <= x * x);
            return positive ? r + x : -r - x;
        }
        Real height = static_cast<Real>(table[layer].first - table[layer - 1].first);
        if (table[layer - 1].first + height * uniform01<Real>(g) < std::exp(-0.5 * x * x))
            return positive ? x : -x;
    }
    throw std::runtime_error("Normal distribution: sampling failed");
}

} // namespace randlib::detail
