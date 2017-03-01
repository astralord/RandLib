#include "UniformRand.h"
#include "../BasicRandGenerator.h"

UniformRand::UniformRand(double minValue, double maxValue) :
    BetaRand(1, 1, minValue, maxValue)
{
}

std::string UniformRand::Name() const
{
    return "Uniform(" + toStringWithPrecision(MinValue()) + ", " + toStringWithPrecision(MaxValue()) + ")";
}

double UniformRand::f(double x) const
{
    return (x < a || x > b) ? 0.0 : bmaInv;
}

double UniformRand::logf(double x) const
{
    return (x < a || x > b) ? -INFINITY : -logBma;
}

double UniformRand::F(double x) const
{
    if (x < a)
        return 0.0;
    return (x > b) ? 1.0 : bmaInv * (x - a);
}

double UniformRand::S(double x) const
{
    if (x < a)
        return 1.0;
    return (x > b) ? 0.0 : bmaInv * (b - x);
}

double UniformRand::Variate() const
{
    return a + StandardVariate() * bma;
}

double UniformRand::Variate(double minValue, double maxValue)
{
    return minValue + StandardVariate() * (maxValue - minValue);
}

double UniformRand::StandardVariate()
{
#ifdef UNIDBLRAND
    /// generates a random number on [0,1) with 53-bit resolution, using 2 32-bit integer variate
    double x;
    unsigned int a, b;
    a = RandGenerator::variate() >> 6; /// Upper 26 bits
    b = RandGenerator::variate() >> 5; /// Upper 27 bits
    x = (a * 134217728.0 + b) / 9007199254740992.0;
    return x;
#elif defined(JLKISS64RAND)
    /// generates a random number on [0,1) with 53-bit resolution, using 64-bit integer variate
    double x;
    unsigned long long a = RandGenerator::variate();
    a = (a >> 12) | 0x3FF0000000000000ULL; /// Take upper 52 bit
    *(reinterpret_cast<unsigned long long *>(&x)) = a; /// Make a double from bits
    return x - 1.0;
#elif defined(UNICLOSEDRAND)
    /// generates a random number on interval [0,1]
    double x = RandGenerator::variate();
    return x / 4294967295.0;
#elif defined(UNIHALFCLOSEDRAND)
    /// generates a random number on interval [0,1)
    double x = RandGenerator::variate();
    return x / 4294967296.0;
#else
    /// generates a random number on interval (0,1)
    double x = RandGenerator::variate();
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

void UniformRand::FitMinimumMLE(const std::vector<double> &sample)
{
    if (!allElementsAreNotBiggerThen(b, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(b)));
    double minVar = *std::min_element(sample.begin(), sample.end());
    SetSupport(minVar, b);
}

void UniformRand::FitMaximumMLE(const std::vector<double> &sample)
{
    if (!allElementsAreNotLessThen(a, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(a)));
    double maxVar = *std::max_element(sample.begin(), sample.end());
    SetSupport(a, maxVar);
}

void UniformRand::FitMLE(const std::vector<double> &sample)
{
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());
    SetSupport(minVar, maxVar);
}

constexpr char UniformRand::TOO_LARGE_A[];
constexpr char UniformRand::TOO_SMALL_B[];

void UniformRand::FitMinimumMM(const std::vector<double> &sample)
{
    if (!allElementsAreNotBiggerThen(b, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(b)));
    double m = sampleMean(sample);
    double leftBound = 2 * m - b;
    if (!allElementsAreNotLessThen(leftBound, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_LARGE_A + toStringWithPrecision(leftBound)));
    SetSupport(leftBound, b);
}

void UniformRand::FitMaximumMM(const std::vector<double> &sample)
{
    if (!allElementsAreNotLessThen(a, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(a)));
    double m = sampleMean(sample);
    double rightBound = 2 * m - a;
    if (!allElementsAreNotBiggerThen(rightBound, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_SMALL_B + toStringWithPrecision(rightBound)));
    SetSupport(a, rightBound);
}

void UniformRand::FitMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double var = sampleVariance(sample, mean);
    double s = std::sqrt(3 * var);
    double leftBound = mean - s, rightBound = mean + s;
    if (!allElementsAreNotLessThen(leftBound, sample))
        throw std::invalid_argument(fitError(WRONG_RETURN, TOO_LARGE_A + toStringWithPrecision(leftBound)));
    if (!allElementsAreNotBiggerThen(rightBound, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_SMALL_B + toStringWithPrecision(rightBound)));
    SetSupport(leftBound, rightBound);
}

void UniformRand::FitMinimumUMVU(const std::vector<double> &sample)
{
    if (!allElementsAreNotBiggerThen(b, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(b)));
    int n = sample.size();
    double minVar = *std::min_element(sample.begin(), sample.end());
    
    /// E[min] = b - n / (n + 1) * (b - a)
    double minVarAdj = (minVar * (n + 1) - b) / n;
    if (!allElementsAreNotLessThen(minVarAdj, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_LARGE_A + toStringWithPrecision(minVarAdj)));
    SetSupport(minVarAdj, b);
}

void UniformRand::FitMaximumUMVU(const std::vector<double> &sample)
{
    if (!allElementsAreNotLessThen(a, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(a)));
    int n = sample.size();
    double maxVar = *std::max_element(sample.begin(), sample.end());

    /// E[max] = (b - a) * n / (n + 1) + a
    double maxVarAdj = (maxVar * (n + 1) - a) / n;
    if (!allElementsAreNotBiggerThen(maxVarAdj, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_SMALL_B + toStringWithPrecision(maxVarAdj)));
    SetSupport(a, maxVarAdj);
}

void UniformRand::FitUMVU(const std::vector<double> &sample)
{
    int n = sample.size();
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());

    /// E[min] = b - n / (n + 1) * (b - a)
    double minVarAdj = (minVar * n - maxVar) / (n - 1);
    /// E[max] = (b - a) * n / (n + 1) + a
    double maxVarAdj = (maxVar * n - minVar) / (n - 1);
    if (!allElementsAreNotLessThen(minVarAdj, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_LARGE_A + toStringWithPrecision(minVarAdj)));
    if (!allElementsAreNotBiggerThen(maxVarAdj, sample))
        throw std::runtime_error(fitError(WRONG_RETURN, TOO_SMALL_B + toStringWithPrecision(maxVarAdj)));
    SetSupport(minVarAdj, maxVarAdj);
}
