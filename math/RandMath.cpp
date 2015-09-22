#include "RandMath.h"

constexpr long double RandMath::factorialTable[];

bool RandMath::areEqual(double a, double b, double eps)
{
    if (a == b)
        return true;
    a = std::fabs(a);
    b = std::fabs(b);
    if (std::fabs(b - a) < eps * std::max(a, b))
        return true;
    return false;
}

long double RandMath::factorialForSmallValue(int n)
{
    int residue = n % 10;
    if (residue <= 5)
    {
        /// go up
        int nPrev = n - residue;
        long double fact = factorialTable[nPrev / 10];
        for (int i = 1; i <= residue; ++i)
            fact *= nPrev + i;
        return fact;
    }

    /// go  down
    int nNext = n - residue + 10;
    double denominator = 1;
    for (int i = 0; i < 10 - residue; ++i)
        denominator *= nNext - i;
    return factorialTable[nNext / 10] / denominator;
}

long double RandMath::factorial(int n)
{
    return (n > maxFactorialTableValue) ? std::tgamma(static_cast<double>(n + 1)) : factorialForSmallValue(n);
}

long double RandMath::doubleFactorial(int n)
{
    long double n_fact = factorial(n);
    if (n & 1) {
        n <<= 1;
        return factorial(n + 1) / (n * n_fact);
    }
    return (1 << n) * n_fact;
}

long double RandMath::binomialCoef(int n, int k)
{
    if (k > n)
        return 0;
    long double n_fact = factorial(n);
    long double k_fact = factorial(k);
    long double n_k_fact = factorial(n - k);
    return n_fact / (k_fact * n_k_fact);
}

long double RandMath::lowerIncGamma(double a, double x)
{
    double sum = 0;
    double term = 1.0 / a;
    int n = 1;
    while (std::fabs(term) > MIN_POSITIVE)
    {
        sum = sum + term;
        term *= x / (a + n);
        ++n;
    }
    return std::pow(x, a) * std::exp(-x) * sum;
}

long double RandMath::upperIncGamma(double a, double x)
{
    double sum = 0;
    double term = 1;
    int n = 1;
    while (std::fabs(term) > MIN_POSITIVE)
    {
        sum = sum + term;
        term *= (a - n) / x;
        ++n;
    }
    return std::pow(x, a - 1) * std::exp(-x) * sum;
}

double RandMath::betaFun(double a, double b)
{
    double sum = a + b;
    if (sum > 30)
    {
        double lgammaA = std::lgamma(a);
        double lgammaB = (a == b) ? lgammaA : std::lgamma(b);
        return std::exp(lgammaA + lgammaB - std::lgamma(sum));
    }

    if (a > b)
        SWAP(a, b);

    double gammaB = std::tgamma(b);
    double res = gammaB / std::tgamma(sum);
    return (a == b) ? res * gammaB : res * std::tgamma(a);
}

double RandMath::regularizedBetaFun(double x, double a, double b)
{
    double upperBoundary = std::min(1.0, std::max(0.0, x));
    if (a != b)
    {
        return integral([a, b] (double t)
        {
            return std::pow(t, a - 1) * std::pow(1 - t, b - 1);
        },
        0, upperBoundary);
    }
    return integral([a, b] (double t)
    {
        return std::pow(t - t * t, a - 1);
    },
    0, upperBoundary);
}

double RandMath::incompleteBetaFun(double x, double a, double b)
{
    // TODO: it doesn't work with a or b == 0!!!
    // look here for effective implementation: http://www.ams.org/journals/mcom/1967-21-100/S0025-5718-1967-0221730-X/S0025-5718-1967-0221730-X.pdf
    return regularizedBetaFun(x, a, b) / betaFun(a, b);
}

long double RandMath::gammaHalf(unsigned k)
{
    if (k & 1)
    {
        unsigned n = (k - 1) >> 1;
        long double res = factorial(k - 1);
        res /= (factorial(n) * (1 << (n << 1)));
        return res * M_SQRTPI;
    }

    return factorial((k >> 1) - 1);
}

long double RandMath::adaptiveSimpsonsAux(const std::function<double (double)> &funPtr, double a, double b,
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

long double RandMath::integral(const std::function<double (double)> &funPtr,
                               double a, double b, double epsilon, int maxRecursionDepth)
{
    if (a > b)
        SWAP(a, b);
    double c = .5 * (a + b), h = (b - a) / 6.0;
    double fa = funPtr(a), fb = funPtr(b), fc = funPtr(c);
    double S = h * (fa + 4 * fc + fb);
    return adaptiveSimpsonsAux(funPtr, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

bool RandMath::findRoot(const std::function<double (double)> &funPtr, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);

    static constexpr int maxIter = 1e5;
    int iter = 0;
    double newRoot = root;
    double step = epsilon + 1;
    double eps = epsilon;
    do {
        root = newRoot;
        double fun0 = funPtr(root);
        double fun1 = funPtr(root + eps);
        if (areEqual(fun1, fun0))
            eps *= 10; /// degeneracy
        else
        {
            if (epsilon < eps)
                eps *= 0.1;
            step = fun0 * eps / (fun1 - fun0);
            newRoot = root - step;
        }
    } while (std::fabs(step) > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// no convergence
        return false;
    return true;
}

bool RandMath::findRoot(const std::function<double (double)> &funPtr, double a, double b, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);

    double fa = funPtr(a);
    if (fa == 0)
    {
        root = a;
        return true;
    }
    double fb = funPtr(b);
    if (fb == 0)
    {
        root = b;
        return true;
    }
    if (fa * fb > 0)
        return false; /// error - the root is not bracketed
    if (std::fabs(fa) < std::fabs(fb))
    {
        SWAP(a, b);
        SWAP(fa, fb);
    }
    double c = a, fc = fa;
    bool mflag = true;
    double s = b, fs = 1, d = 0;
    while (std::fabs(b - a) > epsilon)
    {
        if (!areEqual(fc, fa) && !areEqual(fb, fc))
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
        else
        {
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
        else
            mflag = false;

        fs = funPtr(s);
        if (std::fabs(fs) < epsilon)
        {
            root = s;
            return true;
        }

        d = c;
        c = b;
        fc = fb;

        if (fa * fs < 0)
        {
            b = s;
            fb = fs;
        }
        else
        {
            a = s;
            fa = fs;
        }

        if (std::fabs(fa) < std::fabs(fb))
        {
            SWAP(a, b);
            SWAP(fa, fb);
        }
    }

    root = (std::fabs(fs) < std::fabs(fb)) ? s : b;
    return true;
}

bool RandMath::findMin(const std::function<double (double)> &funPtr, double a, double b, double &root, double epsilon)
{
    if (a > b)
        SWAP(a, b);
    /// golden ratio procedure
    static constexpr double K = 0.5 * (M_SQRT5 - 1);
    double I0 = b - a, I1 = K * I0;
    double xb = a + I1, xa = b - I1;
    double fa = funPtr(xa), fb = funPtr(xb);
    int iter = 0;
    while (++iter < 1e5)
    {
        I1 *= K;
        if (fa >= fb)
        {
            a = xa; xa = xb; xb = a + I1;
            fa = fb; fb = funPtr(xb);
        }
        else
        {
            b = xb; xb = xa; xa = b - I1;
            fb = fa; fa = funPtr(xa);
        }
        if (I1 < epsilon)
        {
            if (fa < fb)
                root = xa;
            else
                root = xb;
            return true;
        }
    }
    return true;

    /// NOT WORKING BRENT
    /*static constexpr double K = 0.5 * (3 - M_SQRT5);
    double x = 0.5 * (a + b), w = x, v = x;
    double fx = funPtr(x), fw = fx, fv = fx;
    double stepCurrent = b - a, stepPrev = stepCurrent;
    double u = a - 1, fu = fx;


    while (std::fabs(a - b) > epsilon)
    {
        double g = stepPrev;
        stepPrev = stepCurrent;

        if (!areEqual(fx, fw) && !areEqual(fx, fv) && !areEqual(fv, fw))
        {
            /// Quadratic approximation
            double diff1 = x - w, diff2 = x - v;
            double fdiff2 = fx - fv, fdiff1 = fx - fw;
            double numerator = diff1 * diff1 * fdiff2;
            numerator -= diff2 * diff2 * fdiff1;
            double denominator = diff1 * fdiff2 - diff2 * fdiff1;
            u = x - 0.5 * numerator / denominator;
        }
        if (a + epsilon <= u && u <= b - epsilon && std::fabs(u - x) < 0.5 * g)
            stepCurrent = std::fabs(u - x); /// accept u
        else
        {
            if (x < 0.5 * (b - a))
            {
                u = x + K * (b - x); /// gold ratio [x, b]
                stepCurrent = b - x;
            }
            else
            {
                u = x - K * (x - a); /// gold ratio [a, x]
                stepCurrent = x - a;
            }
        }

        if (std::fabs(u - x) < epsilon)
        {
            double sign = u - x < 0 ? -1 : 1;
            u = x + sign * epsilon;
        }

        fu = funPtr(u);
        if (fu <= fx)
        {
            if (u >= x)
                a = x;
            else
                b = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        }
        else
        {
            if (u >= x)
                b = u;
            else
                a = u;
            if (fu <= fw || w == x)
            {
                v = w; fv = fw;
                w = u; fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u; fv = fu;
            }
        }
    }

    if (fu < fx)
        root = u;
    else
        root = x;
    return true;*/
}

double RandMath::linearInterpolation(double a, double b, double fa, double fb, double x)
{
    if (b == a)
        return fa;

    double fx = x - a;
    fx /= (b - a);
    fx *= (fb - fa);
    return fx + fa;
}

double RandMath::harmonicNumber(double exponent, int number)
{
    if (exponent < 1 || number < 1)
        return 0;
    double res = 1.0;
    for (int i = 1; i != number; ++i)
        res += std::pow(i + 1, -exponent);
    return res;
}

double RandMath::modifiedBesselFirstKind(double x, int n)
{
    if (n < 0)
        n = -n;
    if (x < 0)
        return 0;

    /// small x
    if (x < 5)
    {
        double halfX = 0.5 * x;
        double halfXSq = halfX * halfX;
        double addon = 1.0;
        double sum = addon;
        double i = 1.0;
        while (std::fabs(addon) > MIN_POSITIVE) {
            addon *= halfXSq;
            addon /= i * i;
            ++i;
            sum += addon;
        };
        return sum * std::pow(halfX, n) / RandMath::factorial(n);
    }

    /// large x
    double addon = 1.0;
    double sum = addon;
    double denominator = 0.125 * x;
    double n2Sq = 4.0 * n * n;
    double i = 1.0, j = i;
    while (std::fabs(addon) > MIN_POSITIVE) {
        double numerator = j * j - n2Sq;
        addon *= numerator * denominator / i;
        sum += addon;
        ++i;
        j += 2;
    };
    return std::exp(x) / std::sqrt(2.0 * M_PI * x) * sum;
}

double RandMath::BernoulliNumber(int n)
{
    std::vector<double> A(n);
    for (int i = 0; i < n; ++i)
    {
        A[i] = 1.0 / (i + 1);
        for (int j = i; j > 1; --j)
        {
            A[j - 1] -= A[j];
            A[j - 1] *= j;
        }
    }
    return A[0];
}

double RandMath::zetaRiemann(int n)
{
    if (n < 0)
    {
        if ((-n) & 1)
            return -BernoulliNumber(n + 1) / (n + 1);
        return 0;
    }

    if (n & 1)
    {
        if (n == 1)
            return INFINITY;
        //TODO:
        return 0;
    }

    double numerator = BernoulliNumber(n);
    numerator *= std::pow(2.0 * M_PI, n);
    double denominator = 2.0 * RandMath::factorial(n);
    double frac = numerator / denominator;
    return ((n >> 1) & 1) ? frac : -frac;
}


