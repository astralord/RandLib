#include "NumericMath.h"

#include <algorithm>
#include <limits>
#include <vector>

namespace RandMath
{

namespace
{

struct Segment
{
    double a, b, estimate, error;
    bool operator<(const Segment &other) const { return error < other.error; }
};

/**
 * @fn gaussKronrod21
 * 21-point Gauss-Kronrod rule with its embedded 10-point Gauss rule (QUADPACK qk21)
 * @param funPtr integrand
 * @param a lower boundary
 * @param b upper boundary
 * @param absEstimate accumulates the integral of |f| over [a, b]
 * @return segment with the Kronrod estimate and |Kronrod - Gauss| as its error
 */
Segment gaussKronrod21(const std::function<double (double)> &funPtr, double a, double b, double &absEstimate)
{
    static constexpr double xgk[11] = {
        0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
        0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
        0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
        0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
        0.294392862701460198131126603103866, 0.148874338981631210884826001129720, 0.0};
    static constexpr double wgk[11] = {
        0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
        0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
        0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
        0.123491976262065851077208980292818, 0.134709217311473325928054001771707,
        0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
        0.149445554002916905664936468389821};
    static constexpr double wg[5] = {
        0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
        0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
        0.295524224714752870173892994651338};
    double center = 0.5 * (a + b), half = 0.5 * (b - a);
    double fc = funPtr(center);
    double kronrod = wgk[10] * fc, gauss = 0.0, absSum = wgk[10] * std::fabs(fc);
    for (int j = 0; j < 10; ++j) {
        double dx = half * xgk[j];
        double f1 = funPtr(center - dx), f2 = funPtr(center + dx);
        kronrod += wgk[j] * (f1 + f2);
        absSum += wgk[j] * (std::fabs(f1) + std::fabs(f2));
        if (j % 2 == 1)
            gauss += wg[j / 2] * (f1 + f2);
    }
    absEstimate += std::fabs(half) * absSum;
    return {a, b, kronrod * half, std::fabs((kronrod - gauss) * half)};
}

}

double integral(const std::function<double (double)> &funPtr, double a, double b, double epsilon, int maxRecursionDepth)
{
    if (a > b)
        return -integral(funPtr, b, a, epsilon, maxRecursionDepth);
    if (a == b)
        return 0.0;
    /// Global adaptive scheme (QUADPACK QAG): always bisect the segment with the largest
    /// error estimate until the total error is below epsilon relative to the integral,
    /// or to a few ulps of the integral of |f| when the integrand cancels
    const int maxSubdivisions = std::max(64, (1 << std::min(maxRecursionDepth, 24)) / 4);
    double absEstimate = 0.0;
    std::vector<Segment> heap{gaussKronrod21(funPtr, a, b, absEstimate)};
    double estimate = heap.front().estimate, error = heap.front().error;
    for (int i = 0; i < maxSubdivisions; ++i) {
        double tolerance = std::max(epsilon * std::fabs(estimate), 50 * std::numeric_limits<double>::epsilon() * absEstimate);
        if (!(error > tolerance) || !std::isfinite(estimate))
            break;
        std::pop_heap(heap.begin(), heap.end());
        Segment worst = heap.back();
        heap.pop_back();
        double middle = 0.5 * (worst.a + worst.b);
        if (!(middle > worst.a && middle < worst.b)) {
            heap.push_back(worst);
            std::push_heap(heap.begin(), heap.end());
            break;
        }
        double newAbs = 0.0;
        Segment left = gaussKronrod21(funPtr, worst.a, middle, newAbs);
        Segment right = gaussKronrod21(funPtr, middle, worst.b, newAbs);
        estimate += left.estimate + right.estimate - worst.estimate;
        error += left.error + right.error - worst.error;
        heap.push_back(left);
        std::push_heap(heap.begin(), heap.end());
        heap.push_back(right);
        std::push_heap(heap.begin(), heap.end());
    }
    /// Re-sum to avoid drift from the running updates
    double sum = 0.0;
    for (const Segment &segment : heap)
        sum += segment.estimate;
    return sum;
}

bool findRootNewtonFirstOrder2d(const std::function<DoublePair (DoublePair)> &funPtr,
                                const std::function<std::tuple<DoublePair, DoublePair> (DoublePair)> &gradPtr,
                                DoublePair &root, long double funTol, long double stepTol)
{
    /// Sanity check
    funTol = funTol > MIN_POSITIVE ? funTol : MIN_POSITIVE;
    stepTol = stepTol > MIN_POSITIVE ? stepTol : MIN_POSITIVE;
    static constexpr int MAX_ITER = 1e5;
    static constexpr double MAX_STEP = 10;
    int iter = 0;
    double step1 = stepTol + 1, step2 = step1;
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
        if (std::min(diffX, relDiffX) < stepTol)
        {
            double diffY1 = fun1- oldFun.first;
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

}
