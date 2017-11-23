#include "CircularDistribution.h"

template < typename RealType >
CircularDistribution<RealType>::CircularDistribution(double location)
{
    SetLocation(location);
}

template < typename RealType >
void CircularDistribution<RealType>::SetLocation(double location)
{
    loc = location;
}
