#pragma once

#include <array>
#include <span>

#include <randlib/continuous/nakagami.hpp>
#include <randlib/continuous/student_t.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/math/beta.hpp>
#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/math/gamma.hpp>
#include <randlib/math/quadrature.hpp>

namespace randlib
{

/**
 * @brief Noncentral Student's t-distribution
 *
 * Notation: X ~ Noncentral-t(ν, δ)
 *
 * If Z is standard normal and V is chi-squared with ν degrees of freedom,
 * X = (Z + δ) / sqrt(V / ν).
 */
template <typename RealType>
class NoncentralTRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double nu = 1;
    double delta = 0;
    double halfNu = 0.5;
    double halfDeltaSquared = 0;
    int seriesCenter = 0;
    StudentTRand<RealType> central{};
    NakagamiRand<RealType> denominator{};

    long double integralPositiveTail(long double x, long double noncentrality, bool upper) const;
    long double integralDensity(long double x) const;
    static long double kummerM(long double a, long double b, long double z);
    long double logPositiveNormalMoment(long double power, long double noncentrality) const;
    long double asymptoticLogTail(long double x, long double noncentrality, bool density) const;
    long double positiveTail(long double x, long double noncentrality, bool upper, bool allowIntegral) const;
    long double probability(long double x, bool upper, bool allowIntegral = true) const;
    long double inverseDenominatorMoment(int order) const;

public:
    explicit NoncentralTRand(double degree = 1, double noncentrality = 0);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double degree, double noncentrality);
    [[nodiscard]] double GetDegree() const noexcept { return nu; }
    [[nodiscard]] double GetNoncentrality() const noexcept { return delta; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/normal.hpp>

namespace randlib
{

template <typename RealType>
NoncentralTRand<RealType>::NoncentralTRand(double degree, double noncentrality)
{
    SetParameters(degree, noncentrality);
}

template <typename RealType>
std::string NoncentralTRand<RealType>::Name() const
{
    return "Noncentral-t(" + this->toStringWithPrecision(nu) + ", " + this->toStringWithPrecision(delta) + ")";
}

template <typename RealType>
void NoncentralTRand<RealType>::SetParameters(double degree, double noncentrality)
{
    if (!(degree > 0) || !std::isfinite(degree))
        throw std::invalid_argument("Noncentral-t distribution: degree should be finite and positive");
    if (!std::isfinite(noncentrality))
        throw std::invalid_argument("Noncentral-t distribution: noncentrality should be finite");
    nu = degree;
    delta = noncentrality;
    halfNu = 0.5 * nu;
    halfDeltaSquared = 0.5 * delta * delta;
    seriesCenter = (std::isfinite(halfDeltaSquared) && halfDeltaSquared < 1000000)
                       ? static_cast<int>(std::floor(halfDeltaSquared))
                       : -1;
    central.SetDegree(nu);
    denominator.SetParameters(halfNu, 1.0);
}

template <typename RealType>
long double NoncentralTRand<RealType>::integralPositiveTail(long double x, long double noncentrality, bool upper) const
{
    const long double shape = halfNu;
    const long double logCoefficient = randlib::consts::ln2 + shape * std::log(shape) - std::lgamma(shape);
    const long double inverseSqrt2 = 1 / std::sqrt(2.0L);
    return math::integrate<long double>(
        [=](long double w) {
            if (!(w > 0))
                return 0.0L;
            const long double normalProbability =
                0.5L * std::erfc((upper ? x * w - noncentrality : noncentrality - x * w) * inverseSqrt2);
            if (normalProbability == 0)
                return 0.0L;
            const long double logDensity = logCoefficient + (2 * shape - 1) * std::log(w) - shape * w * w;
            return normalProbability * std::exp(logDensity);
        },
        0.0L, std::numeric_limits<long double>::infinity(), 2e-13L, 2048);
}

template <typename RealType>
long double NoncentralTRand<RealType>::integralDensity(long double x) const
{
    const long double shape = halfNu;
    const long double noncentrality = delta;
    const long double logCoefficient =
        randlib::consts::ln2 + shape * std::log(shape) - std::lgamma(shape) - 0.5L * std::log(2 * randlib::consts::pi);
    const auto integrand = [=](long double w) {
        if (!(w > 0))
            return 0.0L;
        const long double normalArgument = x * w - noncentrality;
        const long double logIntegrand =
            logCoefficient + 2 * shape * std::log(w) - shape * w * w - 0.5L * normalArgument * normalArgument;
        return std::exp(logIntegrand);
    };
    const long double quadratic = x * x + 2 * shape;
    const long double saddle =
        (x * noncentrality + std::sqrt(x * x * noncentrality * noncentrality + 8 * shape * quadratic)) /
        (2 * quadratic);
    return math::integrate<long double>(integrand, 0.0L, saddle, 2e-13L, 1024) +
           math::integrate<long double>(integrand, saddle, std::numeric_limits<long double>::infinity(), 2e-13L, 1024);
}

template <typename RealType>
long double NoncentralTRand<RealType>::kummerM(long double a, long double b, long double z)
{
    long double term = 1;
    long double sum = 1;
    for (int n = 1; n < 100000; ++n) {
        term *= (a + n - 1) / (b + n - 1) * z / n;
        sum += term;
        if (std::fabs(term) <= 4 * std::numeric_limits<long double>::epsilon() * std::fabs(sum))
            break;
    }
    return sum;
}

template <typename RealType>
long double NoncentralTRand<RealType>::logPositiveNormalMoment(long double power, long double noncentrality) const
{
    const long double saddle = 0.5L * (noncentrality + std::sqrt(noncentrality * noncentrality + 4 * power));
    const long double logPeak = power * std::log(saddle) - 0.5L * (saddle - noncentrality) * (saddle - noncentrality) -
                                0.5L * std::log(2 * randlib::consts::pi);
    const auto integrand = [=](long double value) {
        if (!(value > 0))
            return 0.0L;
        const long double difference = value - noncentrality;
        const long double logValue =
            power * std::log(value) - 0.5L * difference * difference - 0.5L * std::log(2 * randlib::consts::pi);
        return std::exp(logValue - logPeak);
    };
    const long double scaled =
        math::integrate<long double>(integrand, 0.0L, saddle, 2e-13L, 512) +
        math::integrate<long double>(integrand, saddle, std::numeric_limits<long double>::infinity(), 2e-13L, 512);
    return logPeak + std::log(scaled);
}

template <typename RealType>
long double NoncentralTRand<RealType>::asymptoticLogTail(long double x, long double noncentrality, bool density) const
{
    const long double b = halfNu;
    const long double logX = std::log(x);
    std::array<long double, 8> logarithms{};
    long double maximum = -INFINITY;
    int terms = 0;
    for (int k = 0; k < static_cast<int>(logarithms.size()); ++k) {
        const long double power = nu + 2 * k;
        long double logarithm = (b + k) * std::log(b) - std::lgamma(b) - std::lgamma(k + 1.0L) - std::log(b + k) +
                                logPositiveNormalMoment(power, noncentrality) - power * logX;
        if (density)
            logarithm += std::log(power) - logX;
        logarithms[k] = logarithm;
        maximum = std::max(maximum, logarithm);
        terms = k + 1;
        if (k >= 2 && logarithm < maximum - 40)
            break;
    }
    long double scaled = 0;
    for (int k = 0; k < terms; ++k)
        scaled += ((k % 2 == 0) ? 1 : -1) * std::exp(logarithms[k] - maximum);
    return maximum + std::log(scaled);
}

template <typename RealType>
long double NoncentralTRand<RealType>::positiveTail(long double x, long double noncentrality, bool upper,
                                                    bool allowIntegral) const
{
    if (x == 0)
        return 0.5L * std::erfc((upper ? -noncentrality : noncentrality) / std::sqrt(2.0L));
    if (std::isinf(x))
        return upper ? 0.0L : 1.0L;
    const long double lambda = 0.5L * noncentrality * noncentrality;
    if (lambda == 0)
        return upper ? central.S(static_cast<RealType>(x)) : central.F(static_cast<RealType>(x));
    const long double asymptoticCrossover = std::max(100.0L, 10 * std::fabs(noncentrality) + 10);
    if (x > asymptoticCrossover) {
        const long double tail = std::exp(asymptoticLogTail(x, noncentrality, false));
        return upper ? tail : 1 - tail;
    }
    if (seriesCenter < 0)
        return integralPositiveTail(x, noncentrality, upper);

    const long double logRatio = 2 * std::log(x) - std::log(static_cast<long double>(nu));
    const long double logY = -math::softplus(-logRatio);
    const long double log1mY = -math::softplus(logRatio);
    const long double y = std::exp(logY);
    const long double oneMinusY = std::exp(log1mY);
    const long double b = halfNu;
    const int center = seriesCenter;
    long double a1 = center + 0.5L;
    long double a2 = center + 1.0L;
    const long double logBeta1 = math::logBeta(a1, b);
    const long double logBeta2 = math::logBeta(a2, b);
    long double beta1 =
        upper ? math::ibeta(oneMinusY, b, a1, logBeta1, log1mY, logY) : math::ibeta(y, a1, b, logBeta1, logY, log1mY);
    long double beta2 =
        upper ? math::ibeta(oneMinusY, b, a2, logBeta2, log1mY, logY) : math::ibeta(y, a2, b, logBeta2, logY, log1mY);
    long double increment1 = std::exp(a1 * logY + b * log1mY - std::log(a1) - logBeta1);
    long double increment2 = std::exp(a2 * logY + b * log1mY - std::log(a2) - logBeta2);
    const auto centerValue = static_cast<long double>(center);
    const long double logPoisson = (center == 0) ? -lambda
                                                 : -math::stirlingError<long double>(static_cast<std::size_t>(center)) -
                                                       math::bd0<long double>(centerValue, lambda) -
                                                       0.5L * std::log(2 * randlib::consts::pi * centerValue);
    long double poisson = std::exp(logPoisson);
    const long double gammaRatio = std::exp(math::logBeta(centerValue + 1, 0.5L) - 0.5L * randlib::consts::lnpi);
    long double odd = std::copysign(std::sqrt(lambda) * gammaRatio * poisson, noncentrality);
    long double sum = poisson * beta1 + odd * beta2;
    long double compensation = 0;
    long double absoluteSum = std::fabs(poisson * beta1) + std::fabs(odd * beta2);
    const auto addTerm = [&](long double term) {
        const long double corrected = term - compensation;
        const long double next = sum + corrected;
        compensation = (next - sum) - corrected;
        sum = next;
    };

    long double lowerPoisson = poisson;
    long double lowerOdd = odd;
    long double lowerBeta1 = beta1;
    long double lowerBeta2 = beta2;
    long double lowerIncrement1 = increment1;
    long double lowerIncrement2 = increment2;
    long double lowerA1 = a1;
    long double lowerA2 = a2;
    for (int j = center; j > 0; --j) {
        const long double previousIncrement1 = lowerIncrement1 * lowerA1 / (y * (lowerA1 - 1 + b));
        const long double previousIncrement2 = lowerIncrement2 * lowerA2 / (y * (lowerA2 - 1 + b));
        lowerBeta1 += upper ? -previousIncrement1 : previousIncrement1;
        lowerBeta2 += upper ? -previousIncrement2 : previousIncrement2;
        lowerA1 -= 1;
        lowerA2 -= 1;
        lowerIncrement1 = previousIncrement1;
        lowerIncrement2 = previousIncrement2;
        lowerPoisson *= j / lambda;
        lowerOdd *= (j + 0.5L) / lambda;
        const long double term = lowerPoisson * lowerBeta1 + lowerOdd * lowerBeta2;
        addTerm(term);
        absoluteSum += std::fabs(lowerPoisson * lowerBeta1) + std::fabs(lowerOdd * lowerBeta2);
    }

    int decreasing = 0;
    for (int j = center + 1; j < center + 1000000; ++j) {
        poisson *= lambda / j;
        odd *= lambda / (j + 0.5L);
        beta1 += upper ? increment1 : -increment1;
        beta2 += upper ? increment2 : -increment2;
        increment1 *= y * (a1 + b) / (a1 + 1);
        increment2 *= y * (a2 + b) / (a2 + 1);
        a1 += 1;
        a2 += 1;
        const long double term = poisson * beta1 + odd * beta2;
        const long double bound = std::fabs(poisson * beta1) + std::fabs(odd * beta2);
        addTerm(term);
        absoluteSum += bound;
        decreasing = (bound <= 8 * std::numeric_limits<long double>::epsilon() * absoluteSum) ? decreasing + 1 : 0;
        if (decreasing >= 4)
            break;
    }

    long double result = 0.5L * sum;
    if (!upper)
        result += 0.5L * std::erfc(noncentrality / std::sqrt(2.0L));
    const long double integralCrossover = std::max(50.0L, 10 * std::fabs(noncentrality) + 10);
    if (allowIntegral && upper && noncentrality < 0 &&
        (result <= 0 || result <= 256 * std::numeric_limits<long double>::epsilon() * absoluteSum ||
         (result < 1e-8L && x < integralCrossover)))
        result = integralPositiveTail(x, noncentrality, true);
    if (!std::isfinite(result) || result < 0 || result > 1)
        result = integralPositiveTail(x, noncentrality, upper);
    return std::clamp(result, 0.0L, 1.0L);
}

template <typename RealType>
long double NoncentralTRand<RealType>::probability(long double x, bool upper, bool allowIntegral) const
{
    if (x == 0)
        return 0.5L * std::erfc((upper ? -delta : delta) / std::sqrt(2.0L));
    if (x > 0)
        return positiveTail(x, delta, upper, allowIntegral);
    return positiveTail(-x, -delta, !upper, allowIntegral);
}

template <typename RealType>
double NoncentralTRand<RealType>::F(const RealType &x) const
{
    return static_cast<double>(probability(x, false));
}

template <typename RealType>
double NoncentralTRand<RealType>::S(const RealType &x) const
{
    return static_cast<double>(probability(x, true));
}

template <typename RealType>
double NoncentralTRand<RealType>::f(const RealType &x) const
{
    if (!std::isfinite(x))
        return 0;
    if (delta == 0)
        return central.f(x);
    if (x == 0)
        return std::exp(logf(x));
    const long double reflectedNoncentrality = (x > 0) ? delta : -delta;
    const long double asymptoticCrossover = std::max(100.0L, 10 * std::fabs(reflectedNoncentrality) + 10);
    if (std::fabs(static_cast<long double>(x)) > asymptoticCrossover) {
        return static_cast<double>(
            std::exp(asymptoticLogTail(std::fabs(static_cast<long double>(x)), reflectedNoncentrality, true)));
    }
    const long double square = static_cast<long double>(x) * x;
    if (static_cast<long double>(x) * delta < 0) {
        const long double halfDegree = 0.5L * nu;
        const long double z = halfDeltaSquared * square / (square + nu);
        const long double a = 0.5L * (nu + 1);
        const long double summit = 0.5L * (std::sqrt(z * (4 * a + z)) + z);
        if (summit < 40) {
            const long double first = kummerM(a, 0.5L, z);
            long double second = kummerM(halfDegree + 1, 1.5L, z);
            second *= std::exp(std::lgamma(halfDegree + 1) - std::lgamma(a));
            second *= std::sqrt(2.0L) * delta * static_cast<long double>(x) / std::sqrt(square + nu);
            const long double series = first + second;
            if (series > 4 * std::sqrt(std::numeric_limits<long double>::epsilon()) * first) {
                const long double logCoefficient = -halfDeltaSquared - a * std::log1p(square / nu) + std::lgamma(a) -
                                                   std::lgamma(halfDegree) - 0.5L * std::log(nu * randlib::consts::pi);
                return static_cast<double>(std::exp(logCoefficient) * series);
            }
        }
        return static_cast<double>(integralDensity(x));
    }
    const long double scaledX = static_cast<long double>(x) * std::sqrt(1 + 2 / nu);
    const NoncentralTRand<RealType> shifted(nu + 2, delta);
    const long double first = (x > 0) ? probability(x, true, false) : shifted.probability(scaledX, false, false);
    const long double second = (x > 0) ? shifted.probability(scaledX, true, false) : probability(x, false, false);
    const long double difference = first - second;
    if (!(difference > std::sqrt(std::numeric_limits<long double>::epsilon()) * (std::fabs(first) + std::fabs(second))))
        return static_cast<double>(integralDensity(x));
    return static_cast<double>(nu / static_cast<long double>(x) * difference);
}

template <typename RealType>
double NoncentralTRand<RealType>::logf(const RealType &x) const
{
    if (!std::isfinite(x))
        return -INFINITY;
    if (delta == 0)
        return central.logf(x);
    if (x == 0)
        return static_cast<double>(std::lgamma(0.5L * (nu + 1)) - std::lgamma(0.5L * nu) -
                                   0.5L * std::log(static_cast<long double>(nu) * randlib::consts::pi) -
                                   halfDeltaSquared);
    const long double reflectedNoncentrality = (x > 0) ? delta : -delta;
    const long double asymptoticCrossover = std::max(100.0L, 10 * std::fabs(reflectedNoncentrality) + 10);
    if (std::fabs(static_cast<long double>(x)) > asymptoticCrossover) {
        return static_cast<double>(
            asymptoticLogTail(std::fabs(static_cast<long double>(x)), reflectedNoncentrality, true));
    }
    const double density = f(x);
    return (density > 0) ? std::log(density) : -INFINITY;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NoncentralTRand<RealType>::operator()(G &g) const
{
    return (static_cast<RealType>(delta) + detail::standardNormal<RealType>(g)) / denominator(g);
}

template <typename RealType>
void NoncentralTRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    denominator.Sample(g, outputData);
    for (RealType &value : outputData)
        value = (static_cast<RealType>(delta) + detail::standardNormal<RealType>(g)) / value;
}

template <typename RealType>
long double NoncentralTRand<RealType>::inverseDenominatorMoment(int order) const
{
    const long double halfOrder = 0.5L * order;
    const long double halfDegree = halfNu;
    if (order == 2)
        return halfDegree / (halfDegree - 1);
    if (order == 4)
        return halfDegree * halfDegree / ((halfDegree - 1) * (halfDegree - 2));
    const long double logRatio = math::logBeta(halfDegree - halfOrder, halfOrder) - std::lgamma(halfOrder);
    return std::exp(halfOrder * std::log(halfDegree) + logRatio);
}

template <typename RealType>
long double NoncentralTRand<RealType>::Mean() const
{
    return (nu > 1) ? delta * inverseDenominatorMoment(1) : NAN;
}

template <typename RealType>
long double NoncentralTRand<RealType>::Variance() const
{
    if (nu <= 2)
        return (nu > 1) ? INFINITY : NAN;
    const long double first = inverseDenominatorMoment(1);
    const long double second = inverseDenominatorMoment(2);
    const long double dispersion = std::fma(-first, first, second);
    return std::fma(static_cast<long double>(delta) * delta, dispersion, second);
}

template <typename RealType>
long double NoncentralTRand<RealType>::Skewness() const
{
    if (nu <= 3)
        return NAN;
    const long double first = inverseDenominatorMoment(1);
    const long double second = inverseDenominatorMoment(2);
    const long double third = inverseDenominatorMoment(3);
    const long double dispersion = std::fma(-first, first, second);
    const long double thirdDifference = std::fma(-first, second, third);
    const long double thirdCoefficient = std::fma(-2 * first, dispersion, thirdDifference);
    const long double squaredDelta = static_cast<long double>(delta) * delta;
    const long double variance = std::fma(squaredDelta, dispersion, second);
    const long double centralThird = delta * std::fma(squaredDelta, thirdCoefficient, 3 * thirdDifference);
    return centralThird / std::pow(variance, 1.5L);
}

template <typename RealType>
long double NoncentralTRand<RealType>::ExcessKurtosis() const
{
    if (nu <= 4)
        return (nu > 2) ? INFINITY : NAN;
    const long double first = inverseDenominatorMoment(1);
    const long double second = inverseDenominatorMoment(2);
    const long double third = inverseDenominatorMoment(3);
    const long double fourth = inverseDenominatorMoment(4);
    const long double dispersion = std::fma(-first, first, second);
    const long double thirdDifference = std::fma(-first, second, third);
    const long double fourthDifference = std::fma(-first, third, fourth);
    const long double squaredCoefficient = std::fma(-first, thirdDifference, fourthDifference);
    const long double fourthCoefficient =
        std::fma(-3 * first, std::fma(-first, dispersion, thirdDifference), fourthDifference);
    const long double squaredDelta = static_cast<long double>(delta) * delta;
    const long double variance = std::fma(squaredDelta, dispersion, second);
    const long double centralFourth = std::fma(squaredDelta * squaredDelta, fourthCoefficient,
                                               std::fma(6 * squaredDelta, squaredCoefficient, 3 * fourth));
    return centralFourth / (variance * variance) - 3;
}

template <typename RealType>
RealType NoncentralTRand<RealType>::quantileImpl(double p) const
{
    if (delta == 0)
        return central.Quantile(p);
    if (p == static_cast<double>(probability(0, false)))
        return 0;
    const long double mean = Mean();
    const RealType guess = std::isfinite(mean) ? static_cast<RealType>(mean) : static_cast<RealType>(delta);
    return ContinuousDistribution<RealType>::quantileImpl(p, guess);
}

template <typename RealType>
RealType NoncentralTRand<RealType>::quantileImpl1m(double p) const
{
    if (delta == 0)
        return central.Quantile1m(p);
    if (p == static_cast<double>(probability(0, true)))
        return 0;
    const long double mean = Mean();
    const RealType guess = std::isfinite(mean) ? static_cast<RealType>(mean) : static_cast<RealType>(delta);
    return ContinuousDistribution<RealType>::quantileImpl1m(p, guess);
}

template <typename RealType>
std::complex<double> NoncentralTRand<RealType>::CFImpl(double t) const
{
    if (delta == 0)
        return central.CF(t);
    const long double frequency = t;
    const long double noncentrality = delta;
    const long double shape = halfNu;
    const long double logCoefficient = randlib::consts::ln2 + shape * std::log(shape) - std::lgamma(shape);
    const long double linear = 2 * shape - 1;
    const long double squaredMode =
        (linear + std::sqrt(linear * linear + 8 * shape * frequency * frequency)) / (4 * shape);
    const long double mode = std::sqrt(squaredMode);
    const auto envelope = [logCoefficient, linear, shape, frequency](long double w) {
        if (!(w > 0))
            return 0.0L;
        return std::exp(logCoefficient + linear * std::log(w) - shape * w * w - 0.5L * frequency * frequency / (w * w));
    };
    const auto integratePart = [envelope, frequency, noncentrality, mode](bool imaginary) {
        const auto integrand = [envelope, frequency, noncentrality, imaginary](long double w) {
            const long double amplitude = envelope(w);
            if (amplitude == 0)
                return 0.0L;
            const long double phase = frequency * noncentrality / w;
            return amplitude * (imaginary ? std::sin(phase) : std::cos(phase));
        };
        return math::integrate<long double>(integrand, 0.0L, mode, 2e-12L, 1024) +
               math::integrate<long double>(integrand, mode, std::numeric_limits<long double>::infinity(), 2e-12L,
                                            1024);
    };
    return {static_cast<double>(integratePart(false)), static_cast<double>(integratePart(true))};
}

} // namespace randlib
