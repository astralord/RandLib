#include "ExponentialFamily.h"

template< typename T, typename P >
P ExponentialFamily<T, P>::NaturalParameters() const {
    P sourceParameters = SourceParameters();
    return SourceToNatural(sourceParameters);
}

template< typename T, typename P >
double ExponentialFamily<T, P>::ProbabilityMeasure(T x) const
{
    return std::exp(LogProbabilityMeasure(x));
}

template< typename T, typename P >
double ExponentialFamily<T, P>::LogProbabilityMeasure(T x) const
{
    P theta = this->NaturalParameters();
    P t = this->SufficientStatistic(x);
    double y = t * theta;
    y -= this->LogNormalizer(theta);
    y += this->CarrierMeasure(x);
    return y;
}

template< typename T, typename P >
double ExponentialFamily<T, P>::KullbackLeiblerDivergence(P parameters) const
{
    double KL = this->CrossEntropyAdjusted(parameters);
    KL -= this->EntropyAdjusted();
    return KL;
}

template< typename T, typename P >
double ExponentialFamily<T, P>::CrossEntropyAdjusted(P parameters) const
{
    P theta_p = this->NaturalParameters();
    P theta_q = this->SourceToNatural(parameters);
    double H = this->LogNormalizer(theta_q);
    P grad = this->LogNormalizerGradient(theta_p);
    H -= theta_q * grad;
    return H;
}

template< typename T, typename P >
double ExponentialFamily<T, P>::EntropyAdjusted() const
{
    P sourceParameters = this->SourceParameters();
    return CrossEntropyAdjusted(sourceParameters);
}


template class ExponentialFamily<float, double>;
template class ExponentialFamily<double, double>;
template class ExponentialFamily<long double, double>;

template class ExponentialFamily<float, DoublePair>;
template class ExponentialFamily<double, DoublePair>;
template class ExponentialFamily<long double, DoublePair>;

template class ExponentialFamily<int, double>;
template class ExponentialFamily<long int, double>;
template class ExponentialFamily<long long int, double>;

template class ExponentialFamily<int, DoublePair>;
template class ExponentialFamily<long int, DoublePair>;
template class ExponentialFamily<long long int, DoublePair>;

