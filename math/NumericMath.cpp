#include "NumericMath.h"

namespace RandMath
{

/**
 * @fn adaptiveSimpsonsAux
 * auxiliary function for calculation of integral
 * @param funPtr
 * @param a lower boundary
 * @param b upper boundary
 * @param epsilon
 * @param S
 * @param fa
 * @param fb
 * @param fc
 * @param bottom
 * @return
 */
double adaptiveSimpsonsAux(const std::function<double (double)> &funPtr, double a, double b,
                           double epsilon, double S, double fa, double fb, double fc, int bottom)
{
    double c = .5 * (a + b), h = (b - a) / 12.0;
    double d = .5 * (a + c), e = .5 * (c + b);
    double fd = funPtr(d), fe = funPtr(e);
    double Sleft = h * (fa + 4 * fd + fc);
    double Sright = h * (fc + 4 * fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || std::fabs(S2 - S) <= 15.0 * epsilon)
        return S2 + (S2 - S) / 15.0;
    epsilon *= .5;
    --bottom;

    return adaptiveSimpsonsAux(funPtr, a, c, epsilon, Sleft, fa, fc, fd, bottom) +
           adaptiveSimpsonsAux(funPtr, c, b, epsilon, Sright, fc, fb, fe, bottom);
}

double integral(const std::function<double (double)> &funPtr, double a, double b, double epsilon, int maxRecursionDepth)
{
    if (a > b)
        return -integral(funPtr, b, a, epsilon, maxRecursionDepth);
    if (a == b)
        return 0.0;
    double c = .5 * (a + b), h = (b - a) / 6.0;
    double fa = funPtr(a), fb = funPtr(b), fc = funPtr(c);
    double S = h * (fa + 4 * fc + fb);
    return adaptiveSimpsonsAux(funPtr, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

bool findRoot(const std::function<DoublePair (DoublePair)> &funPtr, const std::function<std::tuple<DoublePair, DoublePair> (DoublePair)> &gradPtr,
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
