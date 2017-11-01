#include "UniformRand.h"
#include "../BasicRandGenerator.h"

UniformRand::UniformRand(double minValue, double maxValue) :
    BetaDistribution(1, 1, minValue, maxValue)
{
}

String UniformRand::Name() const
{
    return "Uniform(" + toStringWithPrecision(MinValue()) + ", " + toStringWithPrecision(MaxValue()) + ")";
}

double UniformRand::f(const double & x) const
{
    return (x < a || x > b) ? 0.0 : bmaInv;
}

double UniformRand::logf(const double & x) const
{
    return (x < a || x > b) ? -INFINITY : -logBma;
}

double UniformRand::F(const double & x) const
{
    if (x < a)
        return 0.0;
    return (x > b) ? 1.0 : bmaInv * (x - a);
}

double UniformRand::S(const double & x) const
{
    if (x < a)
        return 1.0;
    return (x > b) ? 0.0 : bmaInv * (b - x);
}

double UniformRand::Variate() const
{
    return a + StandardVariate() * bma;
}

void UniformRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

double UniformRand::Variate(double minValue, double maxValue)
{
    return (minValue < maxValue) ? minValue + StandardVariate() * (maxValue - minValue) : NAN;
}

double UniformRand::StandardVariate()
{
#ifdef UNIDBLRAND
    /// generates a random number on [0,1) with 53-bit resolution, using 2 32-bit integer variate
    double x;
    unsigned int a, b;
    a = RandGenerator::Variate() >> 6; /// Upper 26 bits
    b = RandGenerator::Variate() >> 5; /// Upper 27 bits
    x = (a * 134217728.0 + b) / 9007199254740992.0;
    return x;
#elif defined(JLKISS64RAND)
    /// generates a random number on [0,1) with 53-bit resolution, using 64-bit integer variate
    double x;
    unsigned long long a = RandGenerator::Variate();
    a = (a >> 12) | 0x3FF0000000000000ULL; /// Take upper 52 bit
    *(reinterpret_cast<unsigned long long *>(&x)) = a; /// Make a double from bits
    return x - 1.0;
#elif defined(UNICLOSEDRAND)
    /// generates a random number on interval [0,1]
    double x = RandGenerator::Variate();
    return x / 4294967295.0;
#elif defined(UNIHALFCLOSEDRAND)
    /// generates a random number on interval [0,1)
    double x = RandGenerator::Variate();
    return x / 4294967296.0;
#else
    /// generates a random number on interval (0,1)
    double x = RandGenerator::Variate();
    x += 0.5;
    x /= 4294967296.0;
    return x;
#endif
}

double UniformRand::Mean() const
{
    return 0.5 * (b + a);
}

double UniformRand::Variance() const
{
    return bma * bma / 12;
}

std::complex<double> UniformRand::CFImpl(double t) const
{
    double cosX = std::cos(t * b), sinX = std::sin(t * b);
    double cosY = std::cos(t * a), sinY = std::sin(t * a);
    std::complex<double> numerator(cosX - cosY, sinX - sinY);
    std::complex<double> denominator(0, t * bma);
    return numerator / denominator;
}

double UniformRand::quantileImpl(double p) const
{
    return a + bma * p;
}

double UniformRand::quantileImpl1m(double p) const
{
    return b - bma * p;
}

double UniformRand::Median() const
{
    return 0.5 * (b + a);
}

double UniformRand::Mode() const
{
    /// this can be any value in [a, b]
    return 0.5 * (b + a);
}

double UniformRand::Skewness() const
{
    return 0.0;
}

double UniformRand::ExcessKurtosis() const
{
    return -1.2;
}

double UniformRand::Entropy() const
{
    return (b == a) ? -INFINITY : std::log(bma);
}

constexpr char UniformRand::TOO_LARGE_A[];
constexpr char UniformRand::TOO_SMALL_B[];

void UniformRand::FitMinimum(const std::vector<double> &sample, bool unbiased)
{
    if (!allElementsAreNotBiggerThan(b, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(b)));
    double minVar = *std::min_element(sample.begin(), sample.end());

    if (unbiased == true) {
        int n = sample.size();
        /// E[min] = b - n / (n + 1) * (b - a)
        double minVarAdj = (minVar * (n + 1) - b) / n;
        if (!allElementsAreNotLessThan(minVarAdj, sample))
            throw std::runtime_error(fitErrorDescription(WRONG_RETURN, TOO_LARGE_A + toStringWithPrecision(minVarAdj)));
        SetSupport(minVarAdj, b);
    }
    else {
        SetSupport(minVar, b);
    }
}

void UniformRand::FitMaximum(const std::vector<double> &sample, bool unbiased)
{
    if (!allElementsAreNotLessThan(a, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(a)));
    double maxVar = *std::max_element(sample.begin(), sample.end());

    if (unbiased == true) {
        int n = sample.size();
        /// E[max] = (b - a) * n / (n + 1) + a
        double maxVarAdj = (maxVar * (n + 1) - a) / n;
        if (!allElementsAreNotBiggerThan(maxVarAdj, sample))
            throw std::runtime_error(fitErrorDescription(WRONG_RETURN, TOO_SMALL_B + toStringWithPrecision(maxVarAdj)));
        SetSupport(a, maxVarAdj);
    }
    else {
        SetSupport(a, maxVar);
    }
}

void UniformRand::Fit(const std::vector<double> &sample, bool unbiased)
{
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());
    if (unbiased == true) {
        int n = sample.size();
        /// E[min] = b - n / (n + 1) * (b - a)
        double minVarAdj = (minVar * n - maxVar) / (n - 1);
        /// E[max] = (b - a) * n / (n + 1) + a
        double maxVarAdj = (maxVar * n - minVar) / (n - 1);
        if (!allElementsAreNotLessThan(minVarAdj, sample))
            throw std::runtime_error(fitErrorDescription(WRONG_RETURN, TOO_LARGE_A + toStringWithPrecision(minVarAdj)));
        if (!allElementsAreNotBiggerThan(maxVarAdj, sample))
            throw std::runtime_error(fitErrorDescription(WRONG_RETURN, TOO_SMALL_B + toStringWithPrecision(maxVarAdj)));
        SetSupport(minVarAdj, maxVarAdj);
    }
    else {
        SetSupport(minVar, maxVar);
    }
}
