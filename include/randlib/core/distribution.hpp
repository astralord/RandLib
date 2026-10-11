#pragma once

#include <cstdint>
#include <iomanip>
#include <optional>
#include <random>
#include <span>
#include <sstream>
#include <string>

#include <randlib/math/constants.hpp>
#include <randlib/core/fwd.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/variates.hpp>
#include <randlib/engine/seeding.hpp>
#include <randlib/engine/uniform.hpp>

namespace randlib
{

/**
 * @brief The ProbabilityDistribution class <BR>
 * Abstract class for all probability distributions
 */
template <typename T>
class ProbabilityDistribution
{
    mutable std::optional<DefaultEngine> ownEngine{};

protected:

    /**
     * @brief MAX_ITER_REJECTION
     * upper boundary for maximum amount of iterations in rejection methods
     * one thousand should be enough to be sure there is a bug
     * (or rejection method is too slow to be used)
     */
    static constexpr size_t MAX_ITER_REJECTION = detail::MAX_ITER_REJECTION;

    std::string toStringWithPrecision(const double a_value, const int n = 6) const;

    ProbabilityDistribution() = default;

    /**
     * @fn engine
     * @return the engine of this object, seeded with random_seed() on first use unless Reseed was called
     */
    DefaultEngine &engine() const
    {
        if (!ownEngine)
            ownEngine.emplace(random_seed());
        return *ownEngine;
    }

    /**
     * @fn variate
     * Every concrete distribution implements this by calling its templated operator()
     * @param g engine
     * @return random variable
     */
    virtual T variate(detail::EngineRef g) const = 0;

    /**
     * @fn sample
     * Fills the output with variate(g); distributions with a faster bulk algorithm override it
     * @param g engine
     * @param outputData
     */
    virtual void sample(detail::EngineRef g, std::span<T> outputData) const
    {
        for (T &var : outputData)
            var = variate(g);
    }

public:
    using value_type = T;

    virtual ~ProbabilityDistribution() = default;

    /**
     * @fn Name
     * @return title of distribution, for instance "Normal(0, 1)"
     */
    [[nodiscard]] virtual std::string Name() const = 0;

    /**
     * @fn MinValue
     * @return minimum possible value that can be achieved by random variable
     */
    [[nodiscard]] virtual T MinValue() const = 0;

    /**
     * @fn MaxValue
     * @return maximum possible value that can be achieved by random variable
     */
    [[nodiscard]] virtual T MaxValue() const = 0;

    /**
     * @fn F
     * @param x
     * @return P(X ≤ x)
     */
    [[nodiscard]] virtual double F(const T & x) const = 0;

    /**
     * @fn CumulativeDistributionFunction
     * @param x input vector
     * @param y output vector: y = P(X ≤ x)
     */
    void CumulativeDistributionFunction(std::span<const T> x, std::span<double> y) const;

    /**
     * @fn S
     * @param x
     * @return P(X > x)
     */
    [[nodiscard]] virtual double S(const T & x) const;

    /**
     * @fn SurvivalFunction
     * @param x input vector
     * @param y output vector: y = P(X > x)
     */
    void SurvivalFunction(std::span<const T> x, std::span<double> y) const;

    /**
     * @fn operator()
     * Draws a variate with the given engine through the virtual interface. Concrete distribution
     * classes hide it with a fully templated version, which is faster.
     * @param g any uniform random bit generator, such as std::mt19937 or DefaultEngine
     * @return random variable
     */
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] T operator()(G &g) const { return variate(g); }

    /**
     * @fn Variate
     * Uses the engine of this object, so it is not thread-safe for one object; use operator()
     * with one engine per thread instead
     * @return random variable
     */
    [[nodiscard]] T Variate() const { return variate(engine()); }

    /**
     * @fn Sample
     * Fills the output with variates drawn with the engine of this object
     * @param outputData
     */
    void Sample(std::span<T> outputData) const { sample(engine(), outputData); }

    /**
     * @fn Sample
     * Fills the output with variates drawn with the given engine
     * @param g any uniform random bit generator
     * @param outputData
     */
    template <std::uniform_random_bit_generator G>
    void Sample(G &g, std::span<T> outputData) const { sample(g, outputData); }

    /**
     * @fn Reseed
     * @param seed seed of the engine of this object
     */
    void Reseed(std::uint64_t seed) { ownEngine.emplace(seed); }

protected:
    static constexpr char POSITIVITY_VIOLATION[] = "All elements should be positive";
    static constexpr char NON_NEGATIVITY_VIOLATION[] = "All elements should be non-negative";
    static constexpr char UPPER_LIMIT_VIOLATION[] = "No element should be greater than ";
    static constexpr char LOWER_LIMIT_VIOLATION[] = "No element should be smaller than ";

    /**
     * @fn fitError
     * @param code reason of the failure
     * @param explanation details
     * @return fit_error whose message names this distribution, the reason and the details
     */
    fit_error fitError(FitErrorCode code, const std::string &explanation) const;
};

} // namespace randlib

#include <randlib/math/numeric.hpp>

namespace randlib
{

template < typename T >
std::string ProbabilityDistribution<T>::toStringWithPrecision(const double a_value, const int n) const
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

template < typename T >
void ProbabilityDistribution<T>::CumulativeDistributionFunction(std::span<const T> x, std::span<double> y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = this->F(x[i]);
}

template < typename T >
double ProbabilityDistribution<T>::S(const T &x) const
{
    return 1.0 - this->F(x);
}

template < typename T >
void ProbabilityDistribution<T>::SurvivalFunction(std::span<const T> x, std::span<double> y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = this->S(x[i]);
}

template < typename T >
constexpr char ProbabilityDistribution<T>::POSITIVITY_VIOLATION[];
template < typename T >
constexpr char ProbabilityDistribution<T>::NON_NEGATIVITY_VIOLATION[];
template < typename T >
constexpr char ProbabilityDistribution<T>::UPPER_LIMIT_VIOLATION[];
template < typename T >
constexpr char ProbabilityDistribution<T>::LOWER_LIMIT_VIOLATION[];

template < typename T >
fit_error ProbabilityDistribution<T>::fitError(FitErrorCode code, const std::string &explanation) const
{
    std::string reason;
    switch (code) {
    case FitErrorCode::WrongSample:
        reason = "Sample couldn't be generated by this distribution. ";
        break;
    case FitErrorCode::NotApplicable:
        reason = "Method cannot be applied here. ";
        break;
    case FitErrorCode::WrongReturn:
        reason = "Method returns invalid parameters. ";
        break;
    case FitErrorCode::TooFewElements:
        reason = "Sample is too small. ";
        break;
    case FitErrorCode::WrongLevel:
        reason = "Significance level should be positive and smaller than 1. ";
        break;
    case FitErrorCode::Undefined:
        reason = "Unknown type of error. ";
        break;
    }
    return fit_error(code, this->Name() + ": " + reason + explanation);
}

/// Univariate discrete
/// Univariate continuous
/// Bivariate discrete
/// Bivariate continuous

} // namespace randlib
