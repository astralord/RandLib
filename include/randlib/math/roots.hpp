#pragma once

#include <randlib/math/constants.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace randlib::math
{

namespace detail
{

/// |b - a| no longer resolvable at the relative precision of Real
template <std::floating_point Real>
bool bracketResolved(Real a, Real b)
{
    const Real scale = std::min(std::fabs(a), std::fabs(b));
    return b - a <= 4 * std::numeric_limits<Real>::epsilon() * scale ||
           b - a <= 4 * std::numeric_limits<Real>::denorm_min();
}

template <std::floating_point Real>
Real safeDivide(Real numerator, Real denominator, Real fallback)
{
    if (std::fabs(denominator) < 1 && std::fabs(denominator * std::numeric_limits<Real>::max()) <= std::fabs(numerator))
        return fallback;
    return numerator / denominator;
}

template <std::floating_point Real>
Real secantStep(Real a, Real b, Real fa, Real fb)
{
    const Real tolerance = 5 * std::numeric_limits<Real>::epsilon();
    Real c = a - (fa / (fb - fa)) * (b - a);
    if (c <= a + std::fabs(a) * tolerance || c >= b - std::fabs(b) * tolerance)
        return a + (b - a) / 2;
    return c;
}

template <std::floating_point Real>
Real newtonQuadraticStep(Real a, Real b, Real d, Real fa, Real fb, Real fd, int steps)
{
    Real B = safeDivide(fb - fa, b - a, std::numeric_limits<Real>::max());
    Real A = safeDivide(fd - fb, d - b, std::numeric_limits<Real>::max());
    A = safeDivide(A - B, d - a, Real(0));
    if (A == 0)
        return secantStep(a, b, fa, fb);
    Real c = (std::signbit(A) == std::signbit(fa)) ? a : b;
    for (int i = 0; i < steps; ++i)
        c -= safeDivide(fa + (B + A * (c - b)) * (c - a), B + A * (2 * c - a - b), 1 + c - a);
    if (!(c > a && c < b))
        c = secantStep(a, b, fa, fb);
    return c;
}

template <std::floating_point Real>
Real inverseCubicStep(Real a, Real b, Real d, Real e, Real fa, Real fb, Real fd, Real fe)
{
    Real q11 = (d - e) * fd / (fe - fd);
    Real q21 = (b - d) * fb / (fd - fb);
    Real q31 = (a - b) * fa / (fb - fa);
    Real d21 = (b - d) * fd / (fd - fb);
    Real d31 = (a - b) * fb / (fb - fa);
    Real q22 = (d21 - q11) * fb / (fe - fb);
    Real q32 = (d31 - q21) * fa / (fd - fa);
    Real d32 = (d31 - q21) * fd / (fd - fa);
    Real q33 = (d32 - q22) * fa / (fe - fa);
    Real c = q31 + q32 + q33 + a;
    if (!(c > a && c < b))
        c = newtonQuadraticStep(a, b, d, fa, fb, fd, 3);
    return c;
}

} // namespace detail

/**
 * @fn toms748
 * Bracketing root finder of Alefeld, Potra and Shi (ACM TOMS algorithm 748): inverse cubic and
 * Newton-quadratic interpolation steps, a double-length secant step and bisection whenever the
 * bracket does not shrink by half; it converges like Brent's method at worst and faster on smooth
 * functions
 * @param f continuous function with f(a) and f(b) of opposite signs
 * @param a lower end of the bracket
 * @param b upper end of the bracket
 * @param fa f(a)
 * @param fb f(b)
 * @param maxIterations maximal number of evaluations of f
 * @return the final bracket [a, b] around the root; a == b if f vanishes exactly
 */
template <std::floating_point Real, class F>
std::pair<Real, Real> toms748(const F &f, Real a, Real b, Real fa, Real fb, int maxIterations = 200)
{
    if (a > b) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    if (fa == 0)
        return {a, a};
    if (fb == 0)
        return {b, b};
    if (std::signbit(fa) == std::signbit(fb))
        throw std::invalid_argument("toms748: the root is not bracketed");

    Real d = 0, fd = 0, e = 0, fe = 0;
    int count = maxIterations;
    /// replaces [a, b] by the sub-bracket around c and keeps the dropped end in (d, fd)
    auto bracket = [&f, &a, &b, &fa, &fb, &d, &fd](Real c) {
        const Real tolerance = 2 * std::numeric_limits<Real>::epsilon();
        if (b - a < 2 * tolerance * std::fabs(a))
            c = a + (b - a) / 2;
        else if (c <= a + std::fabs(a) * tolerance)
            c = a + std::fabs(a) * tolerance;
        else if (c >= b - std::fabs(b) * tolerance)
            c = b - std::fabs(b) * tolerance;
        if (!(c > a && c < b))
            c = a + (b - a) / 2;
        Real fc = f(c);
        if (fc == 0) {
            a = b = c;
            fa = fb = 0;
            return;
        }
        if (std::signbit(fa) != std::signbit(fc)) {
            d = b;
            fd = fb;
            b = c;
            fb = fc;
        } else {
            d = a;
            fd = fa;
            a = c;
            fa = fc;
        }
    };
    auto done = [&] { return count <= 0 || fa == 0 || detail::bracketResolved(a, b); };

    bracket(detail::secantStep(a, b, fa, fb));
    --count;
    if (!done()) {
        e = d;
        fe = fd;
        bracket(detail::newtonQuadraticStep(a, b, d, fa, fb, fd, 2));
        --count;
    }
    const Real minDifference = 32 * std::numeric_limits<Real>::min();
    auto distinct = [&] {
        return std::fabs(fa - fb) >= minDifference && std::fabs(fa - fd) >= minDifference &&
               std::fabs(fa - fe) >= minDifference && std::fabs(fb - fd) >= minDifference &&
               std::fabs(fb - fe) >= minDifference && std::fabs(fd - fe) >= minDifference;
    };
    while (!done()) {
        const Real a0 = a, b0 = b;
        Real c = distinct() ? detail::inverseCubicStep(a, b, d, e, fa, fb, fd, fe)
                            : detail::newtonQuadraticStep(a, b, d, fa, fb, fd, 2);
        e = d;
        fe = fd;
        bracket(c);
        --count;
        if (done())
            break;
        c = distinct() ? detail::inverseCubicStep(a, b, d, e, fa, fb, fd, fe)
                       : detail::newtonQuadraticStep(a, b, d, fa, fb, fd, 3);
        bracket(c);
        --count;
        if (done())
            break;
        /// double-length secant step from the end with the smaller residual
        const bool aIsCloser = std::fabs(fa) < std::fabs(fb);
        const Real u = aIsCloser ? a : b, fu = aIsCloser ? fa : fb;
        c = u - 2 * (fu / (fb - fa)) * (b - a);
        if (std::fabs(c - u) > (b - a) / 2)
            c = a + (b - a) / 2;
        e = d;
        fe = fd;
        bracket(c);
        --count;
        if (done())
            break;
        /// bisect when the bracket did not shrink by half
        if (b - a < (b0 - a0) / 2)
            continue;
        e = d;
        fe = fd;
        bracket(a + (b - a) / 2);
        --count;
    }
    return {a, b};
}

/**
 * @fn findRoot
 * @param f continuous function with f(a) and f(b) of opposite signs
 * @param a lower end of the bracket
 * @param b upper end of the bracket
 * @return root of f in [a, b], to full precision (TOMS 748)
 */
template <std::floating_point Real = double, class F>
Real findRoot(const F &f, std::type_identity_t<Real> a, std::type_identity_t<Real> b, int maxIterations = 200)
{
    auto [lo, hi] = toms748<Real>(f, a, b, f(a), f(b), maxIterations);
    return lo + (hi - lo) / 2;
}

/**
 * @fn findRootNewtonFirstOrder
 * Newton's root-finding procedure,
 * using first derivative
 * @param funPtr mapping x |-> (f(x), f'(x))
 * @param root starting point in input and such x that f(x) = 0 in output
 * @param funTol function tolerance
 * @param stepTol step tolerance
 * @return true if success, false otherwise
 */
template <typename RealType>
bool findRootNewtonFirstOrder(const std::function<DoublePair(RealType)> &funPtr, RealType &root,
                              long double funTol = 1e-10, long double stepTol = 1e-6)
{
    /// Sanity check
    funTol = funTol > MIN_POSITIVE ? funTol : MIN_POSITIVE;
    stepTol = stepTol > MIN_POSITIVE ? stepTol : MIN_POSITIVE;
    static constexpr int MAX_ITER = 1e5;
    static constexpr double MAX_STEP = 10;
    int iter = 0;
    double step = stepTol + 1;
    DoublePair y = funPtr(root);
    double fun = y.first;
    double grad = y.second;
    if (std::fabs(fun) < MIN_POSITIVE)
        return true;
    do {
        double alpha = 1.0;
        double oldRoot = root;
        double oldFun = fun;
        step = std::min(MAX_STEP, std::max(-MAX_STEP, fun / grad));
        do {
            root = oldRoot - alpha * step;
            y = funPtr(root);
            fun = y.first;
            grad = y.second;
            if (std::fabs(fun) < MIN_POSITIVE)
                return true;
            alpha *= 0.5;
        } while ((std::fabs(grad) <= MIN_POSITIVE || std::fabs(oldFun) < std::fabs(fun)) && alpha > 0);
        /// Check convergence criteria
        double diffX = std::fabs(root - oldRoot);
        double relDiffX = std::fabs(diffX / oldRoot);
        if (std::min(diffX, relDiffX) < stepTol) {
            double diffY = fun - oldFun;
            double relDiffY = std::fabs(diffY / oldFun);
            if (std::min(std::fabs(fun), relDiffY) < funTol)
                return true;
        }
    } while (++iter < MAX_ITER);
    return false;
}

/**
 * @fn findRootNewtonSecondOrder
 * Newton's root-finding procedure,
 * using first and second derivatives
 * @param funPtr mapping x |-> (f(x), f'(x), f''(x))
 * @param root starting point in input and such x that f(x) = 0 in output
 * @param funTol function tolerance
 * @param stepTol step tolerance
 * @return true if success, false otherwise
 */
template <typename RealType>
bool findRootNewtonSecondOrder(const std::function<DoubleTriplet(RealType)> &funPtr, RealType &root,
                               long double funTol = 1e-10, long double stepTol = 1e-6)
{
    /// Sanity check
    funTol = funTol > MIN_POSITIVE ? funTol : MIN_POSITIVE;
    stepTol = stepTol > MIN_POSITIVE ? stepTol : MIN_POSITIVE;
    static constexpr int MAX_ITER = 1e5;
    static constexpr double MAX_STEP = 10;
    int iter = 0;
    double step = stepTol + 1;
    auto [f, fx, fxx] = funPtr(root);
    if (std::fabs(f) < MIN_POSITIVE)
        return true;
    do {
        double alpha = 1.0;
        double oldRoot = root;
        double oldFun = f;
        double numerator = 2 * f * fx;
        double denominator = 2 * fx * fx - f * fxx;
        step = std::min(MAX_STEP, std::max(-MAX_STEP, numerator / denominator));
        do {
            root = oldRoot - alpha * step;
            std::tie(f, fx, fxx) = funPtr(root);
            if (std::fabs(f) < MIN_POSITIVE)
                return true;
            alpha *= 0.5;
        } while ((std::fabs(fx) <= MIN_POSITIVE || std::fabs(oldFun) < std::fabs(f)) && alpha > 0);
        /// Check convergence criteria
        double diffX = std::fabs(root - oldRoot);
        double relDiffX = std::fabs(diffX / oldRoot);
        if (std::min(diffX, relDiffX) < stepTol) {
            double diffY = f - oldFun;
            double relDiffY = std::fabs(diffY / oldFun);
            if (std::min(std::fabs(f), relDiffY) < funTol)
                return true;
        }
    } while (++iter < MAX_ITER);
    return false;
}

/**
 * @fn findRootNewtonFirstOrder2d
 * Newton's root-finding procedure
 * for 2 functions of 2 parameters,
 * using Jacobian matrix
 * @param funPtr mapping x, y |-> (f(x, y), g(x, y))
 * @param gradPtr x, y |-> (f_x(x, y), f_y(x, y), g_x(x, y)), g_y(x, y))
 * @param root starting point in input and such x that f(x) = 0 in output
 * @param funTol function tolerance
 * @param stepTol step tolerance
 * @return true if success, false otherwise
 */
bool findRootNewtonFirstOrder2d(const std::function<DoublePair(DoublePair)> &funPtr,
                                const std::function<std::tuple<DoublePair, DoublePair>(DoublePair)> &gradPtr,
                                DoublePair &root, long double funTol = 1e-10, long double stepTol = 1e-6);

inline bool findRootNewtonFirstOrder2d(const std::function<DoublePair(DoublePair)> &funPtr,
                                       const std::function<std::tuple<DoublePair, DoublePair>(DoublePair)> &gradPtr,
                                       DoublePair &root, long double funTol, long double stepTol)
{
    /// Sanity check
    funTol = funTol > MIN_POSITIVE ? funTol : MIN_POSITIVE;
    stepTol = stepTol > MIN_POSITIVE ? stepTol : MIN_POSITIVE;
    static constexpr int MAX_ITER = 1e5;
    static constexpr double MAX_STEP = 10;
    int iter = 0;
    double step1;
    double step2;
    DoublePair fun = funPtr(root);
    double fun1 = fun.first;
    double fun2 = fun.second;
    double error = std::max(std::fabs(fun1), std::fabs(fun2));
    if (error < MIN_POSITIVE)
        return true;

    auto [grad1, grad2] = gradPtr(root);
    do {
        double alpha = 1.0;
        DoublePair oldRoot = root;
        DoublePair oldFun = fun;
        double oldError = error;
        double det = grad1.first * grad2.second - grad1.second * grad2.first;
        step1 = std::min(MAX_STEP, std::max(-MAX_STEP, (grad2.second * fun1 - grad1.second * fun2) / det));
        step2 = std::min(MAX_STEP, std::max(-MAX_STEP, (-grad2.first * fun1 + grad1.first * fun2) / det));
        do {
            root.first = oldRoot.first - alpha * step1;
            root.second = oldRoot.second - alpha * step2;
            fun = funPtr(root);
            fun1 = fun.first;
            fun2 = fun.second;
            error = std::max(std::fabs(fun1), std::fabs(fun2));
            if (error < MIN_POSITIVE)
                return true;
            std::tie(grad1, grad2) = gradPtr(root);
            det = grad1.first * grad2.second - grad1.second * grad2.first;
            alpha *= 0.5;
        } while ((std::fabs(det) <= MIN_POSITIVE || oldError < error) && alpha > 0);

        /// Check convergence criteria
        double diffX1 = std::fabs(root.first - oldRoot.first);
        double diffX2 = std::fabs(root.second - oldRoot.second);
        double diffX = std::max(diffX1, diffX2);
        double relDiffX1 = std::fabs(diffX1 / oldRoot.first);
        double relDiffX2 = std::fabs(diffX2 / oldRoot.second);
        double relDiffX = std::max(relDiffX1, relDiffX2);
        if (std::min(diffX, relDiffX) < stepTol) {
            double diffY1 = fun1 - oldFun.first;
            double diffY2 = fun2 - oldFun.second;
            double relDiffY1 = std::fabs(diffY1 / oldFun.first);
            double relDiffY2 = std::fabs(diffY2 / oldFun.second);
            double relDiffY = std::max(relDiffY1, relDiffY2);
            if (std::min(std::max(std::fabs(fun1), std::fabs(fun2)), relDiffY) < funTol)
                return true;
        }
    } while (++iter < MAX_ITER);
    return false;
}

} // namespace randlib::math
