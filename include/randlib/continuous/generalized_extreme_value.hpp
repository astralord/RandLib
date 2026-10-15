#pragma once

#include <randlib/core/continuous.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/math/constants.hpp>

#include <cmath>
#include <limits>
#include <string>

namespace randlib
{

template <typename RealType>
class GeneralizedExtremeValueRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    double shape = 0;
    double scale = 1;
    double location = 0;
    double logScale = 0;

    [[nodiscard]] double standardized(const RealType &x) const { return (x - location) / scale; }
    [[nodiscard]] double supportTerm(const RealType &x) const { return std::fma(shape, standardized(x), 1.0); }
    [[nodiscard]] double supportLog(const RealType &x) const
    {
        const double product = shape * standardized(x);
        return std::fabs(product) < 0.5 ? std::log1p(product) : std::log(supportTerm(x));
    }

public:
    GeneralizedExtremeValueRand(double shapeParameter = 0, double scaleParameter = 1, double locationParameter = 0)
    {
        SetParameters(shapeParameter, scaleParameter, locationParameter);
    }

    void SetParameters(double newShape, double newScale, double newLocation)
    {
        if (!std::isfinite(newShape))
            throw std::invalid_argument("Generalized extreme-value distribution: shape should be finite");
        if (!(newScale > 0) || !std::isfinite(newScale))
            throw std::invalid_argument("Generalized extreme-value distribution: scale should be positive and finite");
        if (!std::isfinite(newLocation))
            throw std::invalid_argument("Generalized extreme-value distribution: location should be finite");
        shape = newShape;
        scale = newScale;
        location = newLocation;
        logScale = std::log(scale);
    }

    [[nodiscard]] double GetShape() const noexcept { return shape; }
    [[nodiscard]] double GetScale() const noexcept { return scale; }
    [[nodiscard]] double GetLocation() const noexcept { return location; }

    [[nodiscard]] std::string Name() const override
    {
        return "GeneralizedExtremeValue(" + this->toStringWithPrecision(shape) + ", " +
               this->toStringWithPrecision(scale) + ", " + this->toStringWithPrecision(location) + ")";
    }

    [[nodiscard]] Support SupportType() const override
    {
        if (shape > 0)
            return Support::RightSemiFinite;
        if (shape < 0)
            return Support::LeftSemiFinite;
        return Support::Infinite;
    }

    [[nodiscard]] RealType MinValue() const override
    {
        return shape > 0 ? static_cast<RealType>(location - scale / shape) : -std::numeric_limits<RealType>::infinity();
    }

    [[nodiscard]] RealType MaxValue() const override
    {
        return shape < 0 ? static_cast<RealType>(location - scale / shape) : std::numeric_limits<RealType>::infinity();
    }

    [[nodiscard]] double logf(const RealType &x) const override
    {
        if (shape == 0) {
            const double z = standardized(x);
            return -logScale - z - std::exp(-z);
        }
        const double support = supportTerm(x);
        if (!(support > 0)) {
            if (support == 0 && shape < 0) {
                if (shape < -1)
                    return INFINITY;
                if (shape == -1)
                    return -logScale;
            }
            return -INFINITY;
        }
        const double logSupport = supportLog(x);
        return -logScale - (1 / shape + 1) * logSupport - std::exp(-logSupport / shape);
    }

    [[nodiscard]] double f(const RealType &x) const override { return std::exp(logf(x)); }

    [[nodiscard]] double F(const RealType &x) const override
    {
        if (shape > 0 && x <= MinValue())
            return 0;
        if (shape < 0 && x >= MaxValue())
            return 1;
        const double exponent = shape == 0 ? std::exp(-standardized(x)) : std::exp(-supportLog(x) / shape);
        return std::exp(-exponent);
    }

    [[nodiscard]] double S(const RealType &x) const override
    {
        if (shape > 0 && x <= MinValue())
            return 1;
        if (shape < 0 && x >= MaxValue())
            return 0;
        const double exponent = shape == 0 ? std::exp(-standardized(x)) : std::exp(-supportLog(x) / shape);
        return -std::expm1(-exponent);
    }

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const
    {
        const RealType exponential = detail::standardExponential<RealType>(g);
        if (shape == 0)
            return static_cast<RealType>(location - scale * std::log(exponential));
        return static_cast<RealType>(location + scale * std::expm1(-shape * std::log(exponential)) / shape);
    }

    [[nodiscard]] long double Mean() const override
    {
        if (shape >= 1)
            return INFINITY;
        if (shape == 0)
            return location + scale * randlib::consts::egamma;
        return location + scale * (std::tgamma(1 - shape) - 1) / shape;
    }

    [[nodiscard]] long double Variance() const override
    {
        if (shape >= 0.5)
            return INFINITY;
        if (shape == 0) {
            const long double value = randlib::consts::pi * scale;
            return value * value / 6;
        }
        const long double first = std::tgamma(1 - shape);
        return scale * scale * (std::tgamma(1 - 2 * shape) - first * first) / (shape * shape);
    }

    [[nodiscard]] RealType Median() const override { return quantileImpl(0.5); }

    [[nodiscard]] RealType Mode() const override
    {
        if (shape <= -1)
            return MaxValue();
        if (shape == 0)
            return static_cast<RealType>(location);
        return static_cast<RealType>(location + scale * std::expm1(-shape * std::log1p(shape)) / shape);
    }

    [[nodiscard]] long double Skewness() const override
    {
        if (shape >= 1.0 / 3)
            return INFINITY;
        if (shape == 0) {
            return 12 * randlib::consts::sqrt2 * randlib::consts::sqrt3 * randlib::consts::apery /
                   (randlib::consts::pi_sq * randlib::consts::pi);
        }
        const long double first = std::tgamma(1 - shape);
        const long double second = std::tgamma(1 - 2 * shape);
        const long double third = std::tgamma(1 - 3 * shape);
        const long double variance = second - first * first;
        const long double centralThird = third - 3 * first * second + 2 * first * first * first;
        return centralThird / std::pow(variance, 1.5L) * (shape > 0 ? 1 : -1);
    }

    [[nodiscard]] long double ExcessKurtosis() const override
    {
        if (shape >= 0.25)
            return INFINITY;
        if (shape == 0)
            return 2.4L;
        const long double first = std::tgamma(1 - shape);
        const long double second = std::tgamma(1 - 2 * shape);
        const long double third = std::tgamma(1 - 3 * shape);
        const long double fourth = std::tgamma(1 - 4 * shape);
        const long double variance = second - first * first;
        const long double centralFourth =
            fourth - 4 * first * third + 6 * first * first * second - 3 * first * first * first * first;
        return centralFourth / (variance * variance) - 3;
    }

    [[nodiscard]] long double Entropy() const { return logScale + randlib::consts::egamma * (1 + shape) + 1; }

private:
    RealType quantileImpl(double p) const override
    {
        const double logExponential = std::log(-std::log(p));
        if (shape == 0)
            return static_cast<RealType>(location - scale * logExponential);
        return static_cast<RealType>(location + scale * std::expm1(-shape * logExponential) / shape);
    }

    RealType quantileImpl1m(double p) const override
    {
        const double logExponential = std::log(-std::log1p(-p));
        if (shape == 0)
            return static_cast<RealType>(location - scale * logExponential);
        return static_cast<RealType>(location + scale * std::expm1(-shape * logExponential) / shape);
    }
};

} // namespace randlib
