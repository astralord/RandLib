#include "RandMath.h"

constexpr long double RandMath::factorialTable[];

bool RandMath::areClose(double a, double b, double eps)
{
    if (a == b)
        return true;
    double fa = std::fabs(a);
    double fb = std::fabs(b);
    return (std::fabs(b - a) < eps * std::max(fa, fb));
}

int RandMath::sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
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

long double RandMath::factorial(double n)
{
    if (n < 0)
        return 0.0;
    return (n > maxFactorialTableValue) ? std::tgamma(n + 1) : factorialForSmallValue(n);
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
    long double n_k_fact;
    if (k == n - k)
        n_k_fact = k_fact;
    else
        n_k_fact = factorial(n - k);
    return n_fact / (k_fact * n_k_fact);
}

double RandMath::digamma(double x)
{
    if (x < 0.0)
        return digamma(1.0 - x) + M_PI / std::tan(M_PI * (1.0 - x));
    double dgam = 0.0;
    if (x > 1000.0)
        return std::log(x) - 0.5 / x;
    while (x > 2.0)
    {
        // TODO: make it faster
        --x;
        dgam += 1.0 / x;
    }
    double y = x - 1.0;
    dgam += y / x - M_EULER;
    // TODO: mininize error by bigger n
    static constexpr int n = 6;
    static constexpr double c[] = {0.64493313, -0.20203181,
                                   0.08209433, -0.03591665,
                                   0.01485925, -0.00472050};
    double r = std::pow(y, n + 1);
    dgam += 0.5 * r;
    for (int i = 0; i != n; ++i)
        dgam += c[i] * (std::pow(y, i + 1) - r);

    if (x < 0.0) /// for x < 0 use Digamma(1-x) = Digamma(x) + pi/tan(pi*x);
        dgam -= M_PI / std::tan(M_PI * x) + 1.0 / x;
    return dgam;
}

double RandMath::trigamma(double x)
{
    if (x < 0.0)
    {
        double z = M_PI / std::sin(M_PI * (1.0 - x));
        return z - digamma(1.0 - x);
    }
    double tgam = 0.0;
    if (x > 1000.0)
    {
        return (x + 0.5) / (x * x);
    }
    while (x > 2.0)
    {
        --x;
        tgam -= 1.0 / (x * x);
    }
    double y = x - 1.0;
    tgam += 1.0 / (x * x);
    // TODO: mininize error by bigger n and don't use same constants twice
    static constexpr int n = 6;
    static constexpr double c[] = {0.64493313, -0.20203181,
                                   0.08209433, -0.03591665,
                                   0.01485925, -0.00472050};
    double r = (n + 1) * std::pow(y, n);
    tgam += 0.5 * r;
    for (int i = 0; i != n; ++i)
        tgam += c[i] * ((i + 1) * std::pow(y, i) - r);
    return tgam;
}

long double RandMath::lowerIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return 0.0;
    if (a == 1)
        return 1.0 - std::exp(-x);
    if (a == 0.5)
        return M_SQRTPI * std::erf(std::sqrt(x));
    if (x > 1.1 && x > a) {
        return std::tgamma(a) - upperIncGamma(a, x);
    }
    return std::exp(logLowerIncGamma(a, x));
}

long double RandMath::logLowerIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return -INFINITY;
    if (a == 1)
        return std::log(1.0 - std::exp(-x));
    if (a == 0.5)
        return 0.5 * M_LNPI + std::log(std::erf(std::sqrt(x)));

    if (x > 1.1 && x > a) {
        double y = std::tgamma(a) - upperIncGamma(a, x);
        return std::log(y);
    }

    long double sum = 0;
    long double term = 1.0 / a;
    double n = a + 1;
    while (std::fabs(term) > MIN_POSITIVE)
    {
        sum = sum + term;
        term *= x / n;
        ++n;
    }
    return a * std::log(x) - x + std::log(sum);
}

long double RandMath::upperIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return std::tgamma(x);
    if (a == 0.5)
        return M_SQRTPI * std::erfc(std::sqrt(x));
    if (a == 1)
        return std::exp(-x);
    if (x <= 1.1 || x <= a) {
        return std::tgamma(a) - lowerIncGamma(a, x);
    }
    return std::exp(logUpperIncGamma(a, x));
}

long double RandMath::logUpperIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return std::lgamma(x);
    if (a == 1)
        return -x;
    if (a == 0.5)
        return 0.5 * M_LNPI + std::log(std::erfc(std::sqrt(x)));
    if (std::max(a, x) < 1e-5)
        return a * std::log(x);
    if (x > 1e4) {
        return (a - 1) * std::log(x) - x;
    }

    if (x <= 1.1 || x <= a) {
        double y = std::tgamma(a) - lowerIncGamma(a, x);
        return std::log(y);
    }

    /// the modified Lentz's method
    double fraction = 1.0 + x - a;
    double C = fraction, D = 0, mult = 1;
    double xmap1 = x - a + 1.0;
    int i = 1;
    do {
        double s = i * (a - i);
        double b = (i << 1) + xmap1;
        D = b + s * D;
        C = b + s / C;
        D = 1.0 / D;
        mult = C * D;
        fraction *= mult;
    } while (std::fabs(mult - 1.0) > MIN_POSITIVE && ++i < 10000);
    return a * std::log(x) - x - std::log(fraction);
}

double RandMath::betaFun(double a, double b)
{
    if (a <= 0 || b <= 0)
        return NAN;
    double lgammaA = std::lgamma(a);
    double lgammaB = (a == b) ? lgammaA : std::lgamma(b);
    return std::exp(lgammaA + lgammaB - std::lgamma(a + b));
}

double RandMath::regularizedBetaFun(double x, double a, double b)
{
    if (a <= 0 || b < 0 || x < 0.0 || x > 1.0)
        return NAN;
    if (b == 0.0 || x == 0.0)
        return 0.0;
    if (x == 1.0)
        return 1.0;
    return incompleteBetaFun(x, a, b) / betaFun(a, b);
}

double RandMath::incompleteBetaFun(double x, double a, double b)
{
    if (a <= 0 || b < 0 || x < 0.0 || x > 1.0) /// if incorrect parameters
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return (b == 0) ? NAN : betaFun(a, b);
    if (a < 1)
    {
        double y = incompleteBetaFun(x, a + 1, b) * (a + b);
        double z = a * std::log(x) + b * std::log(1 - x);
        y += std::exp(z);
        return y / a;
    }
    if (b < 1)
    {
        if (b == 0) /// series expansion
        {
            double denom = a, numen = 1;
            double sum = 1.0 / a, add = 1;
            do {
                numen *= x;
                add = numen / (++denom);
                sum += add;
            } while (add > MIN_POSITIVE);
            return std::pow(x, a) * sum;
        }
        double y = incompleteBetaFun(x, a, b + 1) * (a + b);
        double z = a * std::log(x) + b * std::log(1 - x);
        y -= std::exp(z);
        return y / b;
    }

    double minBound = 0, maxBound = x;
    bool invert = false;
    /// if x > mode
    if (x > (a - 1) / (a + b - 2)) {
        maxBound = 1;
        minBound = x;
        invert = true;
    }
    double y = 0;

    if (a != b)
    {
        y = integral([a, b] (double t)
        {
            double z = (a - 1) * std::log(t) + (b - 1) * std::log(1 - t);
            return std::exp(z);
        },
        minBound, maxBound);
    }
    else {
        y = integral([a, b] (double t)
        {
            return std::pow(t - t * t, a - 1);
        },
        minBound, maxBound);
    }
    return (invert) ? betaFun(a, b) - y : y;
}

long double RandMath::gammaHalf(size_t k)
{
    if (k & 1)
    {
        size_t n = (k - 1) >> 1;
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
        std::swap(a, b);
    if (a == b)
        return 0.0;
    double c = .5 * (a + b), h = (b - a) / 6.0;
    double fa = funPtr(a), fb = funPtr(b), fc = funPtr(c);
    double S = h * (fa + 4 * fc + fb);
    return adaptiveSimpsonsAux(funPtr, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

bool RandMath::findRoot(const std::function<double (double)> &funPtr, const std::function<double (double)> &derPtr, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);
    static constexpr int maxIter = 1e5;
    int iter = 0;
    double step = epsilon + 1;
    double grad = derPtr(root);
    double fun = funPtr(root);
    do {
        double alpha = 1.0;
        double oldRoot = root;
        double oldFun = fun;
        step = fun / grad;
        do {
            root = oldRoot - alpha * step;
            grad = derPtr(root);
            fun = funPtr(root);
            if (std::fabs(fun) < epsilon)
                return true;
            alpha *= 0.5;
        } while ((std::fabs(grad) <= epsilon || std::fabs(oldFun) < std::fabs(fun)) && alpha > 0);
    } while (std::fabs(step) > epsilon && ++iter < maxIter);

    return (iter == maxIter) ? false : true;
}

bool RandMath::findRoot(const std::function<DoublePair (double)> &funPtr, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);
    static constexpr int maxIter = 1e5;
    int iter = 0;
    double step = epsilon + 1;
    DoublePair y = funPtr(root);
    double fun = y.first, grad = y.second;
    do {
        double alpha = 1.0;
        double oldRoot = root;
        double oldFun = fun;
        step = fun / grad;
        do {
            root = oldRoot - alpha * step;
            y = funPtr(root);
            fun = y.first, grad = y.second;
            if (std::fabs(fun) < epsilon)
                return true;
            alpha *= 0.5;
        } while ((std::fabs(grad) <= epsilon || std::fabs(oldFun) < std::fabs(fun)) && alpha > 0);
    } while (std::fabs(step) > epsilon && ++iter < maxIter);

    return (iter == maxIter) ? false : true;
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
        std::swap(a, b);
        std::swap(fa, fb);
    }
    double c = a, fc = fa;
    bool mflag = true;
    double s = b, fs = 1, d = 0;
    while (std::fabs(b - a) > epsilon)
    {
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
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    root = (std::fabs(fs) < std::fabs(fb)) ? s : b;
    return true;
}

bool RandMath::findMin(const std::function<double (double)> &funPtr, double a, double b, double &root, double epsilon)
{
    if (a > b)
        std::swap(a, b);
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
    if (number < 1)
        return 0;
    double res = 1.0;
    for (int i = 2; i <= number; ++i)
        res += std::pow(i, -exponent);
    return res;
}

double RandMath::modifiedBesselFirstKind(double x, double n)
{
    double roundN = std::round(n);
    bool nIsInt = RandMath::areClose(n, roundN);

    if (x < 0) {
        if (nIsInt)
        {
            int nInt = roundN;
            return (nInt % 2) ? -modifiedBesselFirstKind(-x, n) : modifiedBesselFirstKind(-x, n);
        }
        else
            return 0.0;
    }

    if (x == 0) {
        if (n == 0)
            return 1.0;
        return (n > 0 || nIsInt) ? 0.0 : INFINITY;
    }

    if (n == 0.5)
        return std::sqrt(M_2_PI / x) * std::sinh(x);

    if (n == -0.5)
        return std::sqrt(M_2_PI / x) * std::cosh(x);

    /// small x
    if (x < 10)
    {
        double halfX = 0.5 * x;
        double halfXSq = halfX * halfX;
        double addon = 1.0;
        double sum = 1.0;
        double i = 1.0;
        while (std::fabs(addon) > MIN_POSITIVE) {
            addon *= halfXSq;
            addon /= i * (i + n);
            ++i;
            sum += addon;
        }
        double y = std::log(sum);
        y += n * std::log(halfX);
        y -= std::lgamma(n + 1);
        return std::exp(y);
    }

    // large x - divergent sequence!!
    double addon = 1.0;
    double sum = addon;
    double denominator = 0.125 / x;
    double n2Sq = 4.0 * n * n;
    double i = 1.0, j = 1.0;
    while (std::fabs(addon) > MIN_POSITIVE) {
        double numerator = j * j - n2Sq;
        double frac = numerator * denominator / i;
        if (frac > 1)
            break;
        addon *= frac;
        sum += addon;
        ++i;
        j += 2;
    }
    double y = M_PI * x;
    y = std::sqrt(y + y);
    y = std::exp(x) / y;
    return y * sum;
}

double RandMath::modifiedBesselSecondKind(double x, double n)
{
    if (n == 0.5)
        return M_SQRTPI * M_SQRT1_2 * std::exp(-x) / std::sqrt(x);

    double y = modifiedBesselFirstKind(x, -n);
    y -= modifiedBesselFirstKind(x, n);
    y /= std::sin(n * M_PI);
    return 0.5 * M_PI * y;
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

double RandMath::zetaRiemann(double s)
{
    if (s == 1)
        return INFINITY;
    // TODO!! http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf
    if (s < 1)
        return NAN;
    int N = 100;
    double y = harmonicNumber(s, N);
    double NS = std::pow(N, -s);
    return y + N * NS / (s - 1) + 0.5 * NS;
}


