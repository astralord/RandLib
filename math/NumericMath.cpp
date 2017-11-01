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

bool findRoot(const std::function<DoubleTriplet (double)> &funPtr, double &root, double funTol, double stepTol)
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

bool findRoot(const std::function<DoublePair (double)> &funPtr, double &root, double funTol, double stepTol)
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

bool findRoot(const std::function<double (double)> &funPtr, double a, double b, double &root, double epsilon)
{
    /// Sanity check
    epsilon = epsilon > MIN_POSITIVE ? epsilon : MIN_POSITIVE;
    double fa = funPtr(a);
    if (fa == 0) {
        root = a;
        return true;
    }
    double fb = funPtr(b);
    if (fb == 0) {
        root = b;
        return true;
    }
    if (fa * fb > 0) {
        /// error - the root is not bracketed
        return false;
    }
    if (std::fabs(fa) < std::fabs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    double c = a, fc = fa;
    bool mflag = true;
    double s = b, fs = 1, d = 0;
    while (std::fabs(b - a) > epsilon) {
        if (!areClose(fc, fa) && !areClose(fb, fc))
        {
            /// inverse quadratic interpolation
            double numerator = a * fb * fc;
            double denominator = (fa - fb) * (fa - fc);
            s = numerator / denominator;
            numerator = b * fa * fc;
            denominator = (fb - fa) * (fb - fc);
            s += numerator / denominator;
            numerator = c * fa * fb;
            denominator = (fc - fa) * (fc - fb);
            s += numerator / denominator;
        }
        else {
            /// secant method
            s = b - fb * (b - a) / (fb - fa);
        }
        double absDiffSB2 = std::fabs(s - b);
        absDiffSB2 += absDiffSB2;
        double absDiffBC = std::fabs(b - c);
        double absDiffCD = std::fabs(c - d);
        if (s < 0.25 * (3 * a + b) || s > b ||
            (mflag && absDiffSB2 >= absDiffBC) ||
            (!mflag && absDiffSB2 >= absDiffCD) ||
            (mflag && absDiffBC < epsilon) ||
            (!mflag && absDiffCD < epsilon))
        {
            s = 0.5 * (a + b);
            mflag = true;
        }
        else {
            mflag = false;
        }
        fs = funPtr(s);
        if (std::fabs(fs) < epsilon) {
            root = s;
            return true;
        }
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        }
        else {
            a = s;
            fa = fs;
        }
        if (std::fabs(fa) < std::fabs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    root = (std::fabs(fs) < std::fabs(fb)) ? s : b;
    return true;
}

/**
 * @fn parabolicMinimum
 * @param a < b < c
 * @param fa f(a)
 * @param fb f(b)
 * @param fc f(c)
 * @return minimum of interpolated parabola
 */
double parabolicMinimum(double a, double b, double c, double fa, double fb, double fc)
{
    double bma = b - a, cmb = c - b;
    double aux1 = bma * (fb - fc);
    double aux2 = cmb * (fb - fa);
    double numerator = bma * aux1 - cmb * aux2;
    double denominator = aux1 + aux2;
    return b - 0.5 * numerator / denominator;
}

/**
 * @fn findBounds
 * Search of segment that contains the minimum of function
 * @param funPtr mapping x |-> f(x)
 * @param abc such points, that a < b < c, f(a) > f(b) and f(c) > f(b)
 * @param fabc values of a, b and c
 * @param startPoint
 * @return true when segment is found, segment itself and its function values
 */
bool findBounds(const std::function<double (double)> &funPtr, DoubleTriplet &abc, DoubleTriplet &fabc, double startPoint)
{
    static constexpr double K = 0.5 * (M_SQRT5 + 1);
    static constexpr int L = 100;
    double a = startPoint, fa = funPtr(a);
    double b = a + 1.0, fb = funPtr(b);
    double c, fc;
    if (fb < fa) {
        c = b + K * (b - a);
        fc = funPtr(c);
        /// we go to the right
        while (fc < fb) {
            /// parabolic interpolation
            double u = parabolicMinimum(a, b, c, fa, fb, fc);
            double cmb = c - b;
            double fu, uLim = c + L * cmb;
            if (u < c && u > b) {
                fu = funPtr(u);
                if (fu < fc) {
                    abc = std::make_tuple(b, u, c);
                    fabc = std::make_tuple(fb, fu, fc);
                    return true;
                }
                if (fu > fb) {
                    abc = std::make_tuple(a, b, u);
                    fabc = std::make_tuple(fa, fb, fu);
                    return true;
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
        fabc = std::make_tuple(fa, fb, fc);
        return true;
    }
    else {
        c = b; fc = fb;
        b = a; fb = fa;
        a = b - K * (c - b);
        fa = funPtr(a);
        /// go to the left
        while (fa < fb) {
            /// parabolic interpolation
            double u = parabolicMinimum(a, b, c, fa, fb, fc);
            double bma = b - a;
            double fu, uLim = a - L * bma;
            if (u < b && u > a) {
                fu = funPtr(u);
                if (fu < fa) {
                    abc = std::make_tuple(a, u, b);
                    fabc = std::make_tuple(fa, fu, fb);
                    return true;
                }
                if (fu > fb) {
                    abc = std::make_tuple(u, b, c);
                    fabc = std::make_tuple(fu, fb, fc);
                    return true;
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
        fabc = std::make_tuple(fa, fb, fc);
        return true;
    }
}

bool findMin(const std::function<double (double)> &funPtr, const DoubleTriplet & abc, double &fx, double &root, double epsilon)
{
    static constexpr double K = 0.5 * (3 - M_SQRT5);
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
                    u = parabolicMinimum(x, v, w, fx, fv, fw);
                else if (x < w)
                    u = parabolicMinimum(v, x, w, fv, fx, fw);
                else
                    u = parabolicMinimum(v, w, x, fv, fw, fx);
            }
            else {
                if (x < w)
                    u = parabolicMinimum(x, w, v, fx, fv, fw);
                else if (x < v)
                    u = parabolicMinimum(w, x, v, fw, fx, fv);
                else
                    u = parabolicMinimum(w, v, x, fw, fv, fx);
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

bool findMin(const std::function<double (double)> &funPtr, double closePoint, double &root, double epsilon)
{
    DoubleTriplet abc, fabc;
    if (!findBounds(funPtr, abc, fabc, closePoint))
        return false;
    return findMin(funPtr, abc, std::get<1>(fabc), root, epsilon);
}

}
