#include "RandMath.h"

constexpr long double RandMath::factorialTable[];

bool RandMath::areClose(double a, double b, double eps)
{
    if (a == b)
        return true;
    double fa = std::fabs(a);
    double fb = std::fabs(b);
    if (std::fabs(b - a) < eps * std::max(fa, fb))
        return true;
    return false;
}

bool RandMath::sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double RandMath::sum(const QVector<double> &sample)
{
    long double sum = 0.0L;
    for (double var : sample)
        sum += var;
    return sum;
}

double RandMath::sampleMean(const QVector<double> &sample)
{
    int n = sample.size();
    return (n > 0) ? sum(sample) / n : 0;
}

double RandMath::sampleVariance(const QVector<double> &sample, double mean)
{
    int n = sample.size();
    long double deviation = 0.0L;
    for (double var : sample) {
        double diff = (var - mean);
        deviation += diff * diff;
    }
    return (n > 0) ? deviation / n : 0.0;
}

double RandMath::sampleVariance(const QVector<double> &sample)
{
    return sampleVariance(sample, sampleMean(sample));
}

double RandMath::rawMoment(const QVector<double> &sample, int k)
{
    int n = sample.size();
    if (n <= 0 || k < 0)
        return 0.0;
    switch(k) {
        case 0:
            return n;
        case 1:
            return sampleMean(sample);
        default:
        {
            long double sum = 0.0L;
            for (double var : sample)
                sum += std::pow(var, k);
            return sum / n;
        }
    }
}

double RandMath::centralMoment(const QVector<double> &sample, int k, double mean)
{
    int n = sample.size();
    if (n <= 0 || k <= 1)
        return 0.0;
    if (k == 2)
        return sampleVariance(sample, mean);
    long double sum = 0.0L;
    for (double var : sample)
        sum += std::pow(var - mean, k);
    return sum / n;
}

double RandMath::centralMoment(const QVector<double> &sample, int k)
{
    return centralMoment(sample, k, sampleMean(sample));
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
        std::swap(a, b);

    double gammaB = std::tgamma(b);
    double res = gammaB / std::tgamma(sum);
    return (a == b) ? res * gammaB : res * std::tgamma(a);
}

double RandMath::regularizedBetaFun(double x, double a, double b)
{
    if (x < 0.0 || x > 1.0)
        return NAN;
    if (x == 1.0)
        return 1.0;
    if (x == 0.0)
        return 0.0;
    return incompleteBetaFun(x, a, b) / betaFun(a, b);
}

double RandMath::incompleteBetaFun(double x, double a, double b)
{
    if (x < 0.0 || x > 1.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return betaFun(a, b);
    if (a != b)
    {
        return integral([a, b] (double t)
        {
            return std::pow(t, a - 1) * std::pow(1 - t, b - 1);
        },
        0, x);
    }
    return integral([a, b] (double t)
    {
        return std::pow(t - t * t, a - 1);
    },
    0, x);
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
        std::swap(a, b);
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
        step = fun / grad;
        do {
            root = oldRoot - alpha * step;
            grad = derPtr(root);
            fun = funPtr(root);
            if (std::fabs(fun) < epsilon)
                return true;
            alpha *= 0.5;
        } while (grad == 0 && alpha > epsilon);
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
    if (number < 1)
        return 0;
    double res = 1.0;
    for (int i = 2; i <= number; ++i)
        res += std::pow(i, -exponent);
    return res;
}

double RandMath::modifiedBesselFirstKind(double x, int n)
{
    if (n < 0)
        n = -n;

    /// small x
    if (x < 10)
    {
        if (x < 0)
            return 0;
        if (x == 0)
        {
            if (n == 0)
                return 1.0;
            return 0.0;
        }

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
        return sum * std::pow(halfX, n) / RandMath::factorial(n);
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


