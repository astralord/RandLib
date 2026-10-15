#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>

namespace randlib
{

template <typename T, typename P>
class ExponentialFamily
{
public:
    ExponentialFamily() {}
    virtual ~ExponentialFamily() {}

    virtual P SufficientStatistic(T x) const = 0;
    virtual P SourceParameters() const = 0;
    virtual P SourceToNatural(P sourceParameters) const = 0;
    virtual P NaturalParameters() const;

    virtual double LogNormalizer(P theta) const = 0;
    virtual P LogNormalizerGradient(P theta) const = 0;
    virtual double CarrierMeasure(T x) const = 0;

    double ProbabilityMeasure(T x) const;
    double LogProbabilityMeasure(T x) const;

    double KullbackLeiblerDivergence(P parameters) const;
    virtual double CrossEntropyAdjusted(P parameters) const;
    virtual double EntropyAdjusted() const;
};

} // namespace randlib

#include <randlib/core/univariate.hpp>

namespace randlib
{

template <typename T, typename P>
P ExponentialFamily<T, P>::NaturalParameters() const
{
    P sourceParameters = SourceParameters();
    return SourceToNatural(sourceParameters);
}

template <typename T, typename P>
double ExponentialFamily<T, P>::ProbabilityMeasure(T x) const
{
    return std::exp(LogProbabilityMeasure(x));
}

template <typename T, typename P>
double ExponentialFamily<T, P>::LogProbabilityMeasure(T x) const
{
    P theta = this->NaturalParameters();
    P t = this->SufficientStatistic(x);
    double y = t * theta;
    y -= this->LogNormalizer(theta);
    y += this->CarrierMeasure(x);
    return y;
}

template <typename T, typename P>
double ExponentialFamily<T, P>::KullbackLeiblerDivergence(P parameters) const
{
    double KL = this->CrossEntropyAdjusted(parameters);
    KL -= this->EntropyAdjusted();
    return KL;
}

template <typename T, typename P>
double ExponentialFamily<T, P>::CrossEntropyAdjusted(P parameters) const
{
    P theta_p = this->NaturalParameters();
    P theta_q = this->SourceToNatural(parameters);
    double H = this->LogNormalizer(theta_q);
    P grad = this->LogNormalizerGradient(theta_p);
    H -= theta_q * grad;
    return H;
}

template <typename T, typename P>
double ExponentialFamily<T, P>::EntropyAdjusted() const
{
    P sourceParameters = this->SourceParameters();
    return CrossEntropyAdjusted(sourceParameters);
}

} // namespace randlib
