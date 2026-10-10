#include "CircularDistribution.h"

template < typename RealType >
CircularDistribution<RealType>::CircularDistribution(double location)
{
    SetLocation(location);
}

template < typename RealType >
void CircularDistribution<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Circular distribution: location should be finite");
    loc = location;
}

template class CircularDistribution<float>;
template class CircularDistribution<double>;
template class CircularDistribution<long double>;
