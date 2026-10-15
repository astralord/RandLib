#include "stats.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>

namespace randlib_test
{

double kolmogorovSmirnovPValue(double D, size_t n)
{
    double sqrtN = std::sqrt(static_cast<double>(n));
    double lambda = (sqrtN + 0.12 + 0.11 / sqrtN) * D;
    if (lambda < 0.3) {
        /// Q(λ) = 1 - √(2π)/λ Σ exp(-(2k-1)² π² / (8λ²)), accurate for small λ
        const double pi = 3.14159265358979323846;
        double sum = 0.0;
        for (int k = 1; k < 20; ++k) {
            double m = 2 * k - 1;
            sum += std::exp(-m * m * pi * pi / (8 * lambda * lambda));
        }
        return 1.0 - std::sqrt(2 * pi) / lambda * sum;
    }
    double sum = 0.0, sign = 1.0;
    for (int k = 1; k < 200; ++k) {
        double term = std::exp(-2.0 * k * k * lambda * lambda);
        sum += sign * term;
        sign = -sign;
        if (term < 1e-17 * sum)
            break;
    }
    return std::clamp(2.0 * sum, 0.0, 1.0);
}

double kolmogorovSmirnovStatistic(const std::vector<double> &sorted, const std::function<double(double)> &cdf)
{
    /// Compare the empirical distribution function on both sides of every distinct value, using the
    /// left limit of F before it, so that atoms (ties in the sample) are handled correctly
    double n = static_cast<double>(sorted.size()), D = 0.0;
    for (size_t i = 0; i < sorted.size();) {
        size_t j = i;
        while (j < sorted.size() && sorted[j] == sorted[i])
            ++j;
        double v = sorted[i];
        double F = cdf(v), FLeft = cdf(std::nextafter(v, -INFINITY));
        D = std::max({D, std::fabs(FLeft - i / n), std::fabs(F - j / n)});
        i = j;
    }
    return D;
}

namespace
{

double lowerGammaSeries(double a, double x)
{
    double sum = 1.0 / a, term = sum;
    for (int n = 1; n < 100000; ++n) {
        term *= x / (a + n);
        sum += term;
        if (std::fabs(term) < std::fabs(sum) * 1e-17)
            break;
    }
    return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

double upperGammaContinuedFraction(double a, double x)
{
    const double tiny = std::numeric_limits<double>::min() / std::numeric_limits<double>::epsilon();
    double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
    for (int i = 1; i < 100000; ++i) {
        double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < tiny)
            d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny)
            c = tiny;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-16)
            break;
    }
    return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}

struct Segment
{
    double estimate;
    double error;
};

Segment gaussKronrod15(const std::function<double(double)> &f, double a, double b)
{
    static const double xgk[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                                  0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                                  0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                                  0.207784955007898467600689403773245, 0.0};
    static const double wgk[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                                  0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                                  0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                                  0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
    static const double wg[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                                 0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
    double center = 0.5 * (a + b), half = 0.5 * (b - a);
    double fc = f(center);
    double kronrod = wgk[7] * fc, gauss = wg[3] * fc;
    for (int j = 0; j < 7; ++j) {
        double dx = half * xgk[j];
        double sum = f(center - dx) + f(center + dx);
        kronrod += wgk[j] * sum;
        if (j % 2 == 1)
            gauss += wg[j / 2] * sum;
    }
    return {kronrod * half, std::fabs((kronrod - gauss) * half)};
}

double adaptive(const std::function<double(double)> &f, double a, double b, double tolerance)
{
    /// Global adaptive scheme as in QUADPACK's QAG: always bisect the segment with the
    /// largest error estimate, within a fixed budget of subdivisions
    struct Piece
    {
        double a, b;
        Segment s;
    };
    constexpr int MAX_SUBDIVISIONS = 2000;
    std::vector<Piece> pieces{{a, b, gaussKronrod15(f, a, b)}};
    for (int i = 0; i < MAX_SUBDIVISIONS; ++i) {
        double estimate = 0, error = 0;
        size_t worst = 0;
        for (size_t j = 0; j < pieces.size(); ++j) {
            estimate += pieces[j].s.estimate;
            error += pieces[j].s.error;
            if (pieces[j].s.error > pieces[worst].s.error)
                worst = j;
        }
        if (error <= std::max(tolerance * std::fabs(estimate), 1e-300))
            break;
        Piece p = pieces[worst];
        double middle = 0.5 * (p.a + p.b);
        if (!(middle > p.a && middle < p.b))
            break;
        pieces[worst] = {p.a, middle, gaussKronrod15(f, p.a, middle)};
        pieces.push_back({middle, p.b, gaussKronrod15(f, middle, p.b)});
    }
    double estimate = 0;
    for (const Piece &p : pieces)
        estimate += p.s.estimate;
    return estimate;
}

} // namespace

double chiSquaredSurvival(double x, double k)
{
    if (x <= 0.0)
        return 1.0;
    double a = 0.5 * k, y = 0.5 * x;
    return (y < a + 1.0) ? 1.0 - lowerGammaSeries(a, y) : upperGammaContinuedFraction(a, y);
}

double integrate(const std::function<double(double)> &f, double a, double b, double tolerance)
{
    if (a == b)
        return 0.0;
    return adaptive(f, a, b, tolerance);
}

bool relativelyClose(double actual, double expected, double tolerance)
{
    if (std::isnan(actual) || std::isnan(expected))
        return std::isnan(actual) && std::isnan(expected);
    if (actual == expected)
        return true;
    if (std::isinf(actual) || std::isinf(expected))
        return false;
    return std::fabs(actual - expected) <= tolerance * std::max(std::fabs(actual), std::fabs(expected));
}

double relativeError(double actual, double expected)
{
    if (actual == expected)
        return 0.0;
    if (std::isnan(actual) || std::isnan(expected) || std::isinf(actual) || std::isinf(expected))
        return std::numeric_limits<double>::infinity();
    return (expected == 0.0) ? std::fabs(actual) : std::fabs(actual - expected) / std::fabs(expected);
}

void parallelFor(size_t count, const std::function<void(size_t)> &body)
{
    std::atomic<size_t> next{0};
    size_t workers = std::max<size_t>(1, std::min<size_t>(count, std::thread::hardware_concurrency()));
    std::vector<std::thread> threads;
    for (size_t w = 0; w < workers; ++w) {
        threads.emplace_back([&] {
            for (size_t i = next++; i < count; i = next++)
                body(i);
        });
    }
    for (auto &t : threads)
        t.join();
}

} // namespace randlib_test
