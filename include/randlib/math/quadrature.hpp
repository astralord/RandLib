#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace randlib::math
{

namespace detail
{

template <std::floating_point Real>
std::pair<std::complex<Real>, Real> accelerateOscillatorySums(const std::vector<std::complex<Real>> &partialSums)
{
    std::vector<std::complex<Real>> older(partialSums.size() + 1);
    std::vector<std::complex<Real>> previous = partialSums;
    std::complex<Real> best = partialSums.back();
    std::complex<Real> lastEstimate = best;
    Real bestDifference = std::numeric_limits<Real>::infinity();
    const int maxOrder = std::min<int>(20, static_cast<int>(partialSums.size()) - 1);
    for (int order = 1; order <= maxOrder; ++order) {
        std::vector<std::complex<Real>> current(previous.size() - 1);
        for (std::size_t i = 0; i < current.size(); ++i) {
            const std::complex<Real> difference = previous[i + 1] - previous[i];
            const Real scale = std::max({Real{1}, std::abs(previous[i]), std::abs(previous[i + 1])});
            if (std::abs(difference) <= 20 * std::numeric_limits<Real>::epsilon() * scale)
                current[i] = std::complex<Real>(std::numeric_limits<Real>::infinity(), 0);
            else
                current[i] = older[i + 1] + Real{1} / difference;
        }
        if (order % 2 == 0) {
            const std::complex<Real> estimate = current.back();
            if (std::isfinite(estimate.real()) && std::isfinite(estimate.imag())) {
                const Real difference = std::abs(estimate - lastEstimate);
                if (difference < bestDifference) {
                    best = estimate;
                    bestDifference = difference;
                }
                lastEstimate = estimate;
            }
        }
        older = std::move(previous);
        previous = std::move(current);
    }
    return {best, bestDifference};
}

template <std::floating_point Real>
struct Segment
{
    Real a, b, estimate, error;
    bool operator<(const Segment &other) const { return error < other.error; }
};

/**
 * 21-point Gauss-Kronrod rule with its embedded 10-point Gauss rule (QUADPACK qk21)
 * @param f integrand
 * @param a lower boundary
 * @param b upper boundary
 * @param absEstimate accumulates the integral of |f| over [a, b]
 * @return segment with the Kronrod estimate and |Kronrod - Gauss| as its error
 */
template <std::floating_point Real, class F>
Segment<Real> gaussKronrod21(const F &f, Real a, Real b, Real &absEstimate)
{
    static constexpr long double xgk[11] = {0.995657163025808080735527280689003L,
                                            0.973906528517171720077964012084452L,
                                            0.930157491355708226001207180059508L,
                                            0.865063366688984510732096688423493L,
                                            0.780817726586416897063717578345042L,
                                            0.679409568299024406234327365114874L,
                                            0.562757134668604683339000099272694L,
                                            0.433395394129247190799265943165784L,
                                            0.294392862701460198131126603103866L,
                                            0.148874338981631210884826001129720L,
                                            0.0L};
    static constexpr long double wgk[11] = {0.011694638867371874278064396062192L, 0.032558162307964727478818972459390L,
                                            0.054755896574351996031381300244580L, 0.075039674810919952767043140916190L,
                                            0.093125454583697605535065465083366L, 0.109387158802297641899210590325805L,
                                            0.123491976262065851077208980292818L, 0.134709217311473325928054001771707L,
                                            0.142775938577060080797094273138717L, 0.147739104901338491374841515972068L,
                                            0.149445554002916905664936468389821L};
    static constexpr long double wg[5] = {0.066671344308688137593568809893332L, 0.149451349150580593145776339657697L,
                                          0.219086362515982043995534934228163L, 0.269266719309996355091226921569469L,
                                          0.295524224714752870173892994651338L};
    const Real center = (a + b) / 2, half = (b - a) / 2;
    const Real fc = f(center);
    Real kronrod = static_cast<Real>(wgk[10]) * fc, gauss = 0, absSum = static_cast<Real>(wgk[10]) * std::fabs(fc);
    for (int j = 0; j < 10; ++j) {
        const Real dx = half * static_cast<Real>(xgk[j]);
        const Real f1 = f(center - dx), f2 = f(center + dx);
        kronrod += static_cast<Real>(wgk[j]) * (f1 + f2);
        absSum += static_cast<Real>(wgk[j]) * (std::fabs(f1) + std::fabs(f2));
        if (j % 2 == 1)
            gauss += static_cast<Real>(wg[j / 2]) * (f1 + f2);
    }
    absEstimate += std::fabs(half) * absSum;
    return {a, b, kronrod * half, std::fabs((kronrod - gauss) * half)};
}

/// Global adaptive scheme (QUADPACK QAG) on a finite interval
template <std::floating_point Real, class F>
Real integrateFinite(const F &f, Real a, Real b, Real tolerance, int maxSubdivisions)
{
    /// Always bisect the segment with the largest error estimate until the total error is below
    /// the tolerance relative to the integral, or to a few ulps of the integral of |f| when the
    /// integrand cancels
    Real absEstimate = 0;
    std::vector<Segment<Real>> heap{gaussKronrod21(f, a, b, absEstimate)};
    Real estimate = heap.front().estimate, error = heap.front().error;
    /// QUADPACK's roundoff detection: bisections that stop reducing the error mean that rounding
    /// errors of the integrand dominate, and further subdivision cannot reach the tolerance
    int stalledBisections = 0, growingBisections = 0;
    for (int i = 0; i < maxSubdivisions && stalledBisections < 6 && growingBisections < 20; ++i) {
        Real limit = std::max(tolerance * std::fabs(estimate), 50 * std::numeric_limits<Real>::epsilon() * absEstimate);
        if (!(error > limit) || !std::isfinite(estimate))
            break;
        std::pop_heap(heap.begin(), heap.end());
        Segment<Real> worst = heap.back();
        heap.pop_back();
        Real middle = (worst.a + worst.b) / 2;
        if (!(middle > worst.a && middle < worst.b)) {
            heap.push_back(worst);
            std::push_heap(heap.begin(), heap.end());
            break;
        }
        Real newAbs = 0;
        Segment<Real> left = gaussKronrod21(f, worst.a, middle, newAbs);
        Segment<Real> right = gaussKronrod21(f, middle, worst.b, newAbs);
        const Real childEstimate = left.estimate + right.estimate, childError = left.error + right.error;
        if (left.error > 0 && right.error > 0) {
            if (std::fabs(worst.estimate - childEstimate) <= Real(1e-5) * std::fabs(childEstimate) &&
                childError >= Real(0.99) * worst.error)
                ++stalledBisections;
            if (i >= 10 && childError > worst.error)
                ++growingBisections;
        }
        estimate += left.estimate + right.estimate - worst.estimate;
        error += left.error + right.error - worst.error;
        heap.push_back(left);
        std::push_heap(heap.begin(), heap.end());
        heap.push_back(right);
        std::push_heap(heap.begin(), heap.end());
    }
    /// Re-sum to avoid drift from the running updates
    Real sum = 0;
    for (const Segment<Real> &segment : heap)
        sum += segment.estimate;
    return sum;
}

template <std::floating_point Real, class F>
Real evaluateTanhSinh(const F &f, Real x, Real endpointDistance)
{
    if constexpr (std::invocable<F, Real, Real>)
        return static_cast<Real>(f(x, endpointDistance));
    else
        return static_cast<Real>(f(x));
}

template <std::floating_point Real, class F>
std::pair<Real, Real> tanhSinhLevel(const F &f, Real a, Real b, Real step)
{
    const Real halfWidth = (b - a) / 2;
    const Real center = a + halfWidth;
    const Real halfPi = std::numbers::pi_v<Real> / 2;
    const Real centerWeight = halfWidth * halfPi;
    const Real centerValue = evaluateTanhSinh<Real>(f, center, b - center);
    Real sum = centerWeight * centerValue;
    Real absoluteSum = centerWeight * std::fabs(centerValue);

    for (int k = 1; k < 32768; ++k) {
        const Real t = step * static_cast<Real>(k);
        const Real coshT = std::cosh(t);
        const Real transformed = halfPi * std::sinh(t);
        const Real q = std::exp(-2 * transformed);
        if (q == 0)
            break;
        const Real onePlusQ = 1 + q;
        const Real endpointDistance = (b - a) * q / onePlusQ;
        if (endpointDistance < std::numeric_limits<Real>::min())
            break;
        const Real weight = halfWidth * halfPi * coshT * 4 * q / (onePlusQ * onePlusQ);
        if (weight == 0)
            break;

        const Real left = a + endpointDistance;
        const Real right = b - endpointDistance;
        Real leftValue = 0, rightValue = 0;
        if constexpr (std::invocable<F, Real, Real>) {
            leftValue = evaluateTanhSinh<Real>(f, left, -endpointDistance);
            rightValue = evaluateTanhSinh<Real>(f, right, endpointDistance);
        } else {
            if (left != a)
                leftValue = evaluateTanhSinh<Real>(f, left, -endpointDistance);
            if (right != b)
                rightValue = evaluateTanhSinh<Real>(f, right, endpointDistance);
        }
        sum += weight * (leftValue + rightValue);
        absoluteSum += weight * (std::fabs(leftValue) + std::fabs(rightValue));
    }
    return {step * sum, step * absoluteSum};
}

} // namespace detail

/**
 * Global adaptive Gauss-Kronrod quadrature (QUADPACK QAG). Infinite limits are mapped onto (0, 1]
 * by x = a + (1 - t) / t, as in QUADPACK QAGI; the rule never evaluates the endpoints.
 * @param f integrand
 * @param a lower limit, possibly -∞
 * @param b upper limit, possibly +∞
 * @param tolerance relative error tolerance
 * @param maxSubdivisions maximal number of bisections
 * @return ∫_a^b f(x) dx
 */
template <std::floating_point Real = double, class F>
Real integrate(const F &f, std::type_identity_t<Real> a, std::type_identity_t<Real> b,
               std::type_identity_t<Real> tolerance = 1e-11, int maxSubdivisions = 512)
{
    if (a > b)
        return -integrate<Real>(f, b, a, tolerance, maxSubdivisions);
    if (a == b)
        return 0;
    const bool lowerInfinite = std::isinf(a), upperInfinite = std::isinf(b);
    if (!lowerInfinite && !upperInfinite)
        return detail::integrateFinite<Real>(f, a, b, tolerance, maxSubdivisions);
    auto distance = [](Real t) { return (1 - t) / t; };
    if (lowerInfinite && upperInfinite) {
        return detail::integrateFinite<Real>(
            [&f, &distance](Real t) {
                Real x = distance(t);
                return (f(x) + f(-x)) / (t * t);
            },
            0, 1, tolerance, maxSubdivisions);
    }
    if (upperInfinite) {
        return detail::integrateFinite<Real>([&f, &distance, a](Real t) { return f(a + distance(t)) / (t * t); }, 0, 1,
                                             tolerance, maxSubdivisions);
    }
    return detail::integrateFinite<Real>([&f, &distance, b](Real t) { return f(b - distance(t)) / (t * t); }, 0, 1,
                                         tolerance, maxSubdivisions);
}

/**
 * Double-exponential tanh-sinh quadrature on a finite interval. This rule is intended for
 * integrable endpoint singularities. A callable that accepts `(x, endpointDistance)` receives a
 * negative distance `a - x` near the left endpoint and a positive distance `b - x` near the right
 * endpoint, allowing cancellation-free evaluation when `x` rounds to an endpoint.
 * @param f integrand accepting either `f(x)` or `f(x, endpointDistance)`
 * @param a finite lower limit
 * @param b finite upper limit
 * @param tolerance error tolerance relative to the estimated integral of `|f|`
 * @param maxRefinements maximum number of step halvings
 * @return ∫_a^b f(x) dx
 */
template <std::floating_point Real = double, class F>
Real integrateTanhSinh(const F &f, std::type_identity_t<Real> a, std::type_identity_t<Real> b,
                       std::type_identity_t<Real> tolerance = 1e-11, int maxRefinements = 12)
{
    if (!std::isfinite(a) || !std::isfinite(b))
        throw std::invalid_argument("tanh-sinh quadrature requires finite limits");
    if (a > b)
        return -integrateTanhSinh<Real>(f, b, a, tolerance, maxRefinements);
    if (a == b)
        return 0;
    if (!(tolerance > 0) || maxRefinements < 1)
        throw std::invalid_argument("tanh-sinh quadrature requires a positive tolerance and refinement count");

    Real previous = 0;
    for (int level = 0; level <= maxRefinements; ++level) {
        const Real step = std::ldexp(Real{1}, -level);
        const auto [estimate, absoluteEstimate] = detail::tanhSinhLevel<Real>(f, a, b, step);
        const Real roundoff = 50 * std::numeric_limits<Real>::epsilon() * absoluteEstimate;
        const Real target = std::max(tolerance * absoluteEstimate, roundoff);
        if (level >= 2 && std::fabs(estimate - previous) <= target)
            return estimate;
        previous = estimate;
    }
    return previous;
}

namespace detail
{

template <std::floating_point Real, class F, class TailMass>
std::complex<Real> integrateHalfPeriods(const F &amplitude, Real frequency, Real lower, Real tolerance, int maxCycles,
                                        const TailMass &tailMass)
{
    if (!(frequency > 0) || !std::isfinite(frequency) || !std::isfinite(lower) || !(tolerance > 0) || maxCycles < 4)
        throw std::invalid_argument("oscillatory quadrature requires finite limits, positive frequency and tolerance");

    const Real halfPeriod = std::numbers::pi_v<Real> / frequency;
    std::vector<std::complex<Real>> partialSums;
    partialSums.reserve(maxCycles);
    std::complex<Real> sum = 0;
    int negligibleCycles = 0;
    for (int cycle = 0; cycle < maxCycles; ++cycle) {
        const Real a = lower + static_cast<Real>(cycle) * halfPeriod;
        const Real b = a + halfPeriod;
        const Real cosine = integrate<Real>(
            [&amplitude, frequency](Real x) { return static_cast<Real>(amplitude(x)) * std::cos(frequency * x); }, a, b,
            tolerance * 0.1, 128);
        const Real sine = integrate<Real>(
            [&amplitude, frequency](Real x) { return static_cast<Real>(amplitude(x)) * std::sin(frequency * x); }, a, b,
            tolerance * 0.1, 128);
        const std::complex<Real> contribution(cosine, sine);
        sum += contribution;
        partialSums.push_back(sum);

        /// A known bound on the remaining ∫|amplitude| replaces the negligible-cycle heuristic, which
        /// mistakes empty half-periods below the bulk of the amplitude for a converged tail
        const Real remaining = static_cast<Real>(tailMass(b));
        if (std::isnan(remaining)) {
            negligibleCycles = (std::abs(contribution) <= tolerance) ? negligibleCycles + 1 : 0;
            if (negligibleCycles >= 3)
                return sum;
        } else if (remaining <= tolerance) {
            return sum;
        }
        if (!(remaining > Real{0.5}) && partialSums.size() >= 12 && partialSums.size() % 4 == 0) {
            const auto [estimate, error] = accelerateOscillatorySums(partialSums);
            if (error <= tolerance * std::max(Real{1}, std::abs(estimate)))
                return estimate;
        }
    }
    return accelerateOscillatorySums(partialSums).first;
}

} // namespace detail

/**
 * QAWF-style Fourier integral on a semi-infinite interval. The interval is divided into
 * half-periods and the resulting alternating sequence is accelerated with Wynn's epsilon
 * algorithm.
 * @param amplitude non-oscillatory real amplitude
 * @param frequency positive angular frequency
 * @param lower finite lower limit
 * @param tolerance absolute convergence tolerance
 * @param maxCycles maximum number of half-periods
 * @return integral of amplitude(x) exp(i frequency x) from lower to infinity
 */
template <std::floating_point Real = double, class F>
std::complex<Real> integrateOscillatory(const F &amplitude, std::type_identity_t<Real> frequency,
                                        std::type_identity_t<Real> lower = 0,
                                        std::type_identity_t<Real> tolerance = 1e-10, int maxCycles = 512)
{
    return detail::integrateHalfPeriods<Real>(amplitude, frequency, lower, tolerance, maxCycles,
                                              [](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
}

/**
 * QAWF-style Fourier integral of an amplitude with a known tail bound: the half-period walk stops
 * once the remaining mass is below the tolerance, and is accelerated only beyond half of the mass
 * @param amplitude non-oscillatory real amplitude
 * @param frequency positive angular frequency
 * @param lower finite lower limit
 * @param tailMass bound on the integral of |amplitude| from its argument to infinity
 * @param tolerance absolute convergence tolerance
 * @param maxCycles maximum number of half-periods
 * @return integral of amplitude(x) exp(i frequency x) from lower to infinity
 */
template <std::floating_point Real = double, class F, class TailMass>
    requires std::invocable<const TailMass &, Real>
std::complex<Real> integrateOscillatory(const F &amplitude, std::type_identity_t<Real> frequency,
                                        std::type_identity_t<Real> lower, const TailMass &tailMass,
                                        std::type_identity_t<Real> tolerance = 1e-10, int maxCycles = 512)
{
    return detail::integrateHalfPeriods<Real>(amplitude, frequency, lower, tolerance, maxCycles, tailMass);
}

} // namespace randlib::math
