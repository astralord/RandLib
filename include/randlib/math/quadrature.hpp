#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <type_traits>
#include <vector>

namespace randlib::math
{

namespace detail
{

template <std::floating_point Real>
struct Segment
{
    Real a, b, estimate, error;
    bool operator<(const Segment &other) const { return error < other.error; }
};

/**
 * @fn gaussKronrod21
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
    for (int i = 0; i < maxSubdivisions; ++i) {
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

} // namespace detail

/**
 * @fn integrate
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

} // namespace randlib::math
