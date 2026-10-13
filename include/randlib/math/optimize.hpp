#pragma once

#include <randlib/math/constants.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <numeric>
#include <span>
#include <tuple>
#include <vector>

namespace randlib::math
{

/**
 * @fn parabolicMinimum
 * @param a < b < c
 * @param fa f(a)
 * @param fb f(b)
 * @param fc f(c)
 * @return minimum of interpolated parabola
 */
template<typename RealType>
double parabolicMinimum(RealType a, RealType b, RealType c, double fa, double fb, double fc)
{
    RealType bma = b - a, cmb = c - b;
    RealType aux1 = bma * (fb - fc);
    RealType aux2 = cmb * (fb - fa);
    RealType numerator = bma * aux1 - cmb * aux2;
    RealType denominator = aux1 + aux2;
    return b - 0.5 * numerator / denominator;
}

/**
 * @fn findMin
 * Combined Brent's method
 * @param funPtr
 * @param abc lower boundary / middle / upper boundary
 * @param fx funPtr(b)
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return true if success
 */
template<typename RealType>
bool findMin(const std::function<double (RealType)> &funPtr, const Triplet<RealType> & abc, double &fx, RealType &root, double epsilon)
{
    static constexpr double K = 0.5 * (3 - randlib::consts::sqrt5);
    auto [a, x, c] = abc;
    double w = x, v = x, fw = fx, fv = fx;
    double d = c - a, e = d;
    double u = a - 1;
    do {
        double g = e;
        e = d;
        bool acceptParabolicU = false;
        if (x != w && x != v && w != v &&
            fx != fw && fx != fv && fw != fv) {
            if (v < w) {
                if (x < v)
                    u = parabolicMinimum<RealType>(x, v, w, fx, fv, fw);
                else if (x < w)
                    u = parabolicMinimum<RealType>(v, x, w, fv, fx, fw);
                else
                    u = parabolicMinimum<RealType>(v, w, x, fv, fw, fx);
            }
            else {
                if (x < w)
                    u = parabolicMinimum<RealType>(x, w, v, fx, fv, fw);
                else if (x < v)
                    u = parabolicMinimum<RealType>(w, x, v, fw, fx, fv);
                else
                    u = parabolicMinimum<RealType>(w, v, x, fw, fv, fx);
            }
            double absumx = std::fabs(u - x);
            if (u >= a + epsilon && u <= c - epsilon && absumx < 0.5 * g) {
                acceptParabolicU = true; /// accept u
                d = absumx;
            }
        }

        if (!acceptParabolicU) {
            /// use golden ratio instead of parabolic approximation
            if (x < 0.5 * (c + a)) {
                d = c - x;
                u = x + K * d; /// golden ratio [x, c]
            }
            else {
                d = x - a;
                u = x - K * d; /// golden ratio [a, x]
            }
        }

        if (std::fabs(u - x) < epsilon) {
            u = x + epsilon * sign(u - x); /// setting the closest distance between u and x
        }

        double fu = funPtr(u);
        if (fu <= fx) {
            if (u >= x)
                a = x;
            else
                c = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        }
        else {
            if (u >= x)
                c = u;
            else
                a = u;
            if (fu <= fw || w == x) {
                v = w; w = u;
                fv = fw; fw = fu;
            }
            else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    } while (0.49 * (c - a) > epsilon);
    root = x;
    return true;
}

/**
 * @fn findMin
 * Combined Brent's method
 * @param funPtr
 * @param closePoint point that is nearby minimum
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return true if success
 */
template<typename RealType>
bool findMin(const std::function<double (RealType)> &funPtr, RealType closePoint, RealType &root, long double epsilon = 1e-8)
{
    Triplet<RealType> abc;
    static constexpr double K = 0.5 * (randlib::consts::sqrt5 + 1);
    static constexpr int L = 100;
    double a = closePoint, fa = funPtr(a);
    double b = a + 1.0, fb = funPtr(b);
    double c, fc;
    if (fb < fa) {
        c = b + K * (b - a);
        fc = funPtr(c);
        /// we go to the right
        while (fc < fb) {
            /// parabolic interpolation
            double u = parabolicMinimum<RealType>(a, b, c, fa, fb, fc);
            double cmb = c - b;
            double fu, uLim = c + L * cmb;
            if (u < c && u > b) {
                fu = funPtr(u);
                if (fu < fc) {
                    abc = std::make_tuple(b, u, c);
                    return findMin(funPtr, abc, fu, root, epsilon);
                }
                if (fu > fb) {
                    abc = std::make_tuple(a, b, u);
                    return findMin(funPtr, abc, fb, root, epsilon);
                }
                u = c + K * cmb;
                fu = funPtr(u);
            }
            else if (u > c && u < uLim) {
                fu = funPtr(u);
                if (fu < fc) {
                    b = c; c = u; u = c + K * cmb;
                    fb = fc, fc = fu, fu = funPtr(u);
                }
            }
            else if (u > uLim) {
                u = uLim;
                fu = funPtr(u);
            }
            else {
                u = c + K * cmb;
                fu = funPtr(u);
            }
            a = b; b = c; c = u;
            fa = fb; fb = fc; fc = fu;
        }
        abc = std::make_tuple(a, b, c);
        return findMin(funPtr, abc, fb, root, epsilon);
    }
    else {
        c = b; fc = fb;
        b = a; fb = fa;
        a = b - K * (c - b);
        fa = funPtr(a);
        /// go to the left
        while (fa < fb) {
            /// parabolic interpolation
            double u = parabolicMinimum<RealType>(a, b, c, fa, fb, fc);
            double bma = b - a;
            double fu, uLim = a - L * bma;
            if (u < b && u > a) {
                fu = funPtr(u);
                if (fu < fa) {
                    abc = std::make_tuple(a, u, b);
                    return findMin(funPtr, abc, fu, root, epsilon);
                }
                if (fu > fb) {
                    abc = std::make_tuple(u, b, c);
                    return findMin(funPtr, abc, fb, root, epsilon);
                }
                u = a - K * bma;
                fu = funPtr(u);
            }
            else if (u < a && u > uLim) {
                fu = funPtr(u);
                if (fu < fa) {
                    b = a; a = u; u = a - K * bma;
                    fb = fa, fa = fu, fu = funPtr(u);
                }
            }
            else if (u < uLim) {
                u = uLim;
                fu = funPtr(u);
            }
            else {
                u = a - K * bma;
                fu = funPtr(u);
            }
            c = b; b = a; a = u;
            fc = fb; fb = fa; fa = fu;
        }
        abc = std::make_tuple(a, b, c);
        return findMin(funPtr, abc, fb, root, epsilon);
    }
}

/**
 * @brief The MinimizeResult struct
 * Outcome of a multivariate minimization
 */
template <std::floating_point Real = double>
struct MinimizeResult
{
    std::vector<Real> x; ///< point of the minimum
    Real value = 0; ///< function value there
    int iterations = 0; ///< number of iterations
    bool converged = false; ///< false if the iteration limit was reached first
};

/**
 * @fn nelderMead
 * Downhill simplex method of Nelder and Mead with the dimension-dependent coefficients of Gao and Han
 * (2012). The initial simplex steps 5% along each coordinate (0.00025 for zero coordinates), as in
 * MATLAB's fminsearch.
 * @param f function of a std::span<const Real>
 * @param start starting point
 * @param tolerance the iteration stops when the function values and the coordinates of the simplex
 * agree to this relative tolerance
 * @param maxIterations iteration limit
 */
template <std::floating_point Real = double, class F>
MinimizeResult<Real> nelderMead(const F &f, std::vector<Real> start, Real tolerance = 1e-10, int maxIterations = 20000)
{
    const std::size_t n = start.size();
    const Real dimension = static_cast<Real>(std::max<std::size_t>(n, 1));
    const Real reflection = 1, expansion = 1 + 2 / dimension, contraction = Real(0.75) - 1 / (2 * dimension),
               shrinkage = 1 - 1 / dimension;
    std::vector<std::vector<Real>> simplex(n + 1, start);
    for (std::size_t i = 0; i < n; ++i)
        simplex[i + 1][i] = (start[i] != 0) ? start[i] * Real(1.05) : Real(0.00025);
    std::vector<Real> values(n + 1);
    for (std::size_t i = 0; i <= n; ++i)
        values[i] = f(std::span<const Real>(simplex[i]));

    std::vector<std::size_t> order(n + 1);
    auto point = [n] (const std::vector<Real> &from, const std::vector<Real> &to, Real t) {
        std::vector<Real> x(n);
        for (std::size_t j = 0; j < n; ++j)
            x[j] = from[j] + t * (to[j] - from[j]);
        return x;
    };
    MinimizeResult<Real> result;
    for (result.iterations = 0; result.iterations < maxIterations; ++result.iterations) {
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&values] (std::size_t i, std::size_t j) { return values[i] < values[j]; });
        const std::size_t best = order.front(), worst = order.back(), secondWorst = order[n > 0 ? n - 1 : 0];

        Real valueSpread = 0, pointSpread = 0, scale = 0;
        for (std::size_t i = 0; i <= n; ++i) {
            valueSpread = std::max(valueSpread, std::fabs(values[i] - values[best]));
            for (std::size_t j = 0; j < n; ++j) {
                pointSpread = std::max(pointSpread, std::fabs(simplex[i][j] - simplex[best][j]));
                scale = std::max(scale, std::fabs(simplex[best][j]));
            }
        }
        if (valueSpread <= tolerance * (std::fabs(values[best]) + tolerance) && pointSpread <= tolerance * (scale + 1)) {
            result.converged = true;
            break;
        }

        std::vector<Real> centroid(n, 0);
        for (std::size_t i : order) {
            if (i == worst)
                continue;
            for (std::size_t j = 0; j < n; ++j)
                centroid[j] += simplex[i][j] / dimension;
        }
        std::vector<Real> reflected = point(centroid, simplex[worst], -reflection);
        Real reflectedValue = f(std::span<const Real>(reflected));
        if (reflectedValue < values[best]) {
            std::vector<Real> expanded = point(centroid, simplex[worst], -reflection * expansion);
            Real expandedValue = f(std::span<const Real>(expanded));
            if (expandedValue < reflectedValue) {
                simplex[worst] = expanded;
                values[worst] = expandedValue;
            }
            else {
                simplex[worst] = reflected;
                values[worst] = reflectedValue;
            }
            continue;
        }
        if (reflectedValue < values[secondWorst]) {
            simplex[worst] = reflected;
            values[worst] = reflectedValue;
            continue;
        }
        const bool outside = reflectedValue < values[worst];
        std::vector<Real> contracted = outside ? point(centroid, reflected, contraction)
                                               : point(centroid, simplex[worst], contraction);
        Real contractedValue = f(std::span<const Real>(contracted));
        if (contractedValue < (outside ? reflectedValue : values[worst])) {
            simplex[worst] = contracted;
            values[worst] = contractedValue;
            continue;
        }
        for (std::size_t i = 0; i <= n; ++i) {
            if (i == best)
                continue;
            simplex[i] = point(simplex[best], simplex[i], shrinkage);
            values[i] = f(std::span<const Real>(simplex[i]));
        }
    }
    std::size_t best = std::min_element(values.begin(), values.end()) - values.begin();
    result.x = simplex[best];
    result.value = values[best];
    return result;
}

/**
 * @fn numericalGradient
 * @param f function of a std::span<const Real>
 * @param x point
 * @return central-difference gradient with steps proportional to the cube root of the machine epsilon
 */
template <std::floating_point Real = double, class F>
std::vector<Real> numericalGradient(const F &f, const std::vector<Real> &x)
{
    const Real h0 = std::cbrt(std::numeric_limits<Real>::epsilon());
    std::vector<Real> gradient(x.size()), y = x;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const Real h = h0 * std::max<Real>(std::fabs(x[i]), 1);
        y[i] = x[i] + h;
        const Real up = f(std::span<const Real>(y));
        y[i] = x[i] - h;
        const Real down = f(std::span<const Real>(y));
        y[i] = x[i];
        gradient[i] = (up - down) / (2 * h);
    }
    return gradient;
}

/**
 * @fn bfgs
 * Quasi-Newton minimization of Broyden, Fletcher, Goldfarb and Shanno with central-difference
 * gradients and a backtracking line search under the Armijo condition
 * @param f smooth function of a std::span<const Real>
 * @param start starting point
 * @param tolerance the iteration stops when the largest gradient component is below this tolerance
 * times max(1, |f|), or when a step no longer changes f
 * @param maxIterations iteration limit
 */
template <std::floating_point Real = double, class F>
MinimizeResult<Real> bfgs(const F &f, std::vector<Real> start, Real tolerance = 1e-8, int maxIterations = 1000)
{
    const std::size_t n = start.size();
    auto identity = [n] {
        std::vector<std::vector<Real>> m(n, std::vector<Real>(n, 0));
        for (std::size_t i = 0; i < n; ++i)
            m[i][i] = 1;
        return m;
    };
    auto dot = [] (const std::vector<Real> &u, const std::vector<Real> &v) {
        return std::inner_product(u.begin(), u.end(), v.begin(), Real(0));
    };
    std::vector<std::vector<Real>> H = identity();
    std::vector<Real> x = std::move(start);
    Real value = f(std::span<const Real>(x));
    std::vector<Real> gradient = numericalGradient<Real>(f, x);
    MinimizeResult<Real> result;
    for (result.iterations = 0; result.iterations < maxIterations; ++result.iterations) {
        Real largest = 0;
        for (Real g : gradient)
            largest = std::max(largest, std::fabs(g));
        if (largest <= tolerance * std::max<Real>(1, std::fabs(value))) {
            result.converged = true;
            break;
        }
        std::vector<Real> direction(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                direction[i] -= H[i][j] * gradient[j];
        Real slope = dot(gradient, direction);
        if (!(slope < 0)) {
            H = identity();
            for (std::size_t i = 0; i < n; ++i)
                direction[i] = -gradient[i];
            slope = -dot(gradient, gradient);
        }
        Real step = 1, newValue = value;
        std::vector<Real> next(n);
        for (int k = 0; k < 60; ++k, step /= 2) {
            for (std::size_t i = 0; i < n; ++i)
                next[i] = x[i] + step * direction[i];
            newValue = f(std::span<const Real>(next));
            if (newValue <= value + Real(1e-4) * step * slope)
                break;
        }
        if (!(newValue < value)) {
            result.converged = newValue == value;
            break;
        }
        std::vector<Real> newGradient = numericalGradient<Real>(f, next);
        std::vector<Real> s(n), y(n);
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = next[i] - x[i];
            y[i] = newGradient[i] - gradient[i];
        }
        const Real sy = dot(s, y);
        if (sy > 0) {
            /// H <- (I - ρ s yᵀ) H (I - ρ y sᵀ) + ρ s sᵀ
            const Real rho = 1 / sy;
            std::vector<Real> Hy(n, 0);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    Hy[i] += H[i][j] * y[j];
            const Real yHy = dot(y, Hy);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    H[i][j] += rho * ((1 + rho * yHy) * s[i] * s[j] - Hy[i] * s[j] - s[i] * Hy[j]);
        }
        x = std::move(next);
        value = newValue;
        gradient = std::move(newGradient);
    }
    result.x = std::move(x);
    result.value = value;
    return result;
}

} // namespace randlib::math
