#include "ExponentialRand.h"
#include "UniformRand.h"
#include "../BasicRandGenerator.h"

String ExponentialRand::Name() const
{
    return "Exponential(" + toStringWithPrecision(GetRate()) + ")";
}

double ExponentialRand::f(const double & x) const
{
    return (x < 0.0) ? 0.0 : beta * std::exp(-beta * x);
}

double ExponentialRand::logf(const double & x) const
{
    return (x < 0.0) ? -INFINITY : logBeta - beta * x;
}

double ExponentialRand::F(const double & x) const
{
    return (x > 0.0) ? -std::expm1(-beta * x) : 0.0;
}

double ExponentialRand::S(const double & x) const
{
    return (x > 0.0) ? std::exp(-beta * x) : 1.0;
}

double ExponentialRand::Variate() const
{
    return theta * StandardVariate(localRandGenerator);
}

void ExponentialRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

double ExponentialRand::StandardVariate(RandGenerator &randGenerator)
{
    /// Ziggurat algorithm
    int iter = 0;
    do {
        int stairId = randGenerator.Variate() & 255;
        /// Get horizontal coordinate
        double x = UniformRand::StandardVariate(randGenerator) * ziggurat[stairId].second;
        if (x < ziggurat[stairId + 1].second) /// if we are under the upper stair - accept
            return x;
        if (stairId == 0) /// if we catch the tail
            return ziggurat[1].second + StandardVariate(randGenerator);
        long double height = ziggurat[stairId].first - ziggurat[stairId - 1].first;
        if (ziggurat[stairId - 1].first + height * UniformRand::StandardVariate(randGenerator) < std::exp(-x)) /// if we are under the curve - accept
            return x;
        /// rejection - go back
    } while (++iter <= MAX_ITER_REJECTION);
    /// fail due to some error
    return NAN;
}

std::complex<double> ExponentialRand::CFImpl(double t) const
{
    return 1.0 / std::complex<double>(1.0, -theta * t);
}

double ExponentialRand::Entropy() const
{
    return 1.0 - logBeta;
}

double ExponentialRand::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return std::exp(RandMath::lfact(n) - n * logBeta);
}
