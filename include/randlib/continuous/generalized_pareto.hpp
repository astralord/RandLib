#pragma once

#include <randlib/core/continuous.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/exponential.hpp>

#include <cmath>
#include <limits>
#include <string>

namespace randlib
{

template <typename RealType>
class GeneralizedParetoRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    double shape = 0;
    double scale = 1;
    double location = 0;
    double logScale = 0;

    [[nodiscard]] double supportLog(double z) const
    {
        const double product = shape * z;
        return std::fabs(product) < 0.5 ? std::log1p(product) : std::log(std::fma(shape, z, 1.0));
    }

public:
    GeneralizedParetoRand(double shapeParameter = 0, double scaleParameter = 1, double locationParameter = 0)
    {
        SetParameters(shapeParameter, scaleParameter, locationParameter);
    }

    void SetParameters(double newShape, double newScale, double newLocation)
    {
        if (!std::isfinite(newShape))
            throw std::invalid_argument("Generalized Pareto distribution: shape should be finite");
        if (!(newScale > 0) || !std::isfinite(newScale))
            throw std::invalid_argument("Generalized Pareto distribution: scale should be positive and finite");
        if (!std::isfinite(newLocation))
            throw std::invalid_argument("Generalized Pareto distribution: location should be finite");
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
        return "GeneralizedPareto(" + this->toStringWithPrecision(shape) + ", " + this->toStringWithPrecision(scale) +
               ", " + this->toStringWithPrecision(location) + ")";
    }

    [[nodiscard]] Support SupportType() const override
    {
        return shape < 0 ? Support::Finite : Support::RightSemiFinite;
    }
    [[nodiscard]] RealType MinValue() const override { return static_cast<RealType>(location); }
    [[nodiscard]] RealType MaxValue() const override
    {
        return shape < 0 ? static_cast<RealType>(location - scale / shape) : std::numeric_limits<RealType>::infinity();
    }

    [[nodiscard]] double logf(const RealType &x) const override
    {
        if (x < location || (shape < 0 && x > MaxValue()))
            return -INFINITY;
        const double z = (x - location) / scale;
        if (shape == 0)
            return -logScale - z;
        const double support = std::fma(shape, z, 1.0);
        if (!(support > 0)) {
            if (support == 0) {
                if (shape < -1)
                    return INFINITY;
                if (shape == -1)
                    return -logScale;
            }
            return -INFINITY;
        }
        return -logScale - (1 / shape + 1) * supportLog(z);
    }

    [[nodiscard]] double f(const RealType &x) const override { return std::exp(logf(x)); }

    [[nodiscard]] double S(const RealType &x) const override
    {
        if (x <= location)
            return 1;
        if (shape < 0 && x >= MaxValue())
            return 0;
        const double z = (x - location) / scale;
        return shape == 0 ? std::exp(-z) : std::exp(-supportLog(z) / shape);
    }

    [[nodiscard]] double F(const RealType &x) const override
    {
        if (x <= location)
            return 0;
        if (shape < 0 && x >= MaxValue())
            return 1;
        const double z = (x - location) / scale;
        const double logSurvival = shape == 0 ? -z : -supportLog(z) / shape;
        return -std::expm1(logSurvival);
    }

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const
    {
        const RealType exponential = detail::standardExponential<RealType>(g);
        if (shape == 0)
            return static_cast<RealType>(location + scale * exponential);
        return static_cast<RealType>(location + scale * std::expm1(shape * exponential) / shape);
    }

    [[nodiscard]] long double Mean() const override { return shape < 1 ? location + scale / (1 - shape) : INFINITY; }

    [[nodiscard]] long double Variance() const override
    {
        if (shape >= 0.5)
            return INFINITY;
        const long double denominator = (1 - shape) * (1 - shape) * (1 - 2 * shape);
        return scale * scale / denominator;
    }

    [[nodiscard]] RealType Median() const override { return quantileImpl(0.5); }

    [[nodiscard]] RealType Mode() const override { return shape < -1 ? MaxValue() : static_cast<RealType>(location); }

    [[nodiscard]] long double Skewness() const override
    {
        if (shape >= 1.0 / 3)
            return INFINITY;
        return 2 * (1 + shape) * std::sqrt(1 - 2 * shape) / (1 - 3 * shape);
    }

    [[nodiscard]] long double ExcessKurtosis() const override
    {
        if (shape >= 0.25)
            return INFINITY;
        const long double numerator = 3 * (1 - 2 * shape) * (2 * shape * shape + shape + 3);
        return numerator / ((1 - 3 * shape) * (1 - 4 * shape)) - 3;
    }

    [[nodiscard]] long double Entropy() const { return logScale + shape + 1; }

private:
    RealType quantileImpl(double p) const override
    {
        const double exponential = -std::log1p(-p);
        if (shape == 0)
            return static_cast<RealType>(location + scale * exponential);
        return static_cast<RealType>(location + scale * std::expm1(shape * exponential) / shape);
    }

    RealType quantileImpl1m(double p) const override
    {
        const double exponential = -std::log(p);
        if (shape == 0)
            return static_cast<RealType>(location + scale * exponential);
        return static_cast<RealType>(location + scale * std::expm1(shape * exponential) / shape);
    }
};

} // namespace randlib
