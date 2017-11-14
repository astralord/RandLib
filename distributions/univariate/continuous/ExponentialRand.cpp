#include "ExponentialRand.h"
#include "UniformRand.h"
#include "../BasicRandGenerator.h"

long double ExponentialRand::stairWidth[257] = {0};
long double ExponentialRand::stairHeight[256] = {0};
bool ExponentialRand::dummy = ExponentialRand::SetupTables();

String ExponentialRand::Name() const
{
    return "Exponential(" + toStringWithPrecision(GetRate()) + ")";
}

bool ExponentialRand::SetupTables()
{
    /// Set up ziggurat tables
    static constexpr long double A = 3.9496598225815571993e-3l; /// area under rectangle
    /// coordinates of the implicit rectangle in base layer
    stairHeight[0] = 0.00045413435384149675l; /// exp(-x1);
    stairWidth[0] = 8.697117470131049720307l; /// A / stairHeight[0];
    /// implicit value for the top layer
    stairWidth[256] = 0;
    stairWidth[1] = x1;
    stairHeight[1] = 0.0009672692823271745203l;
    for (size_t i = 2; i < 256; ++i) {
        /// such y_i that f(x_{i+1}) = y_i
        stairWidth[i] = -std::log(stairHeight[i - 1]);
        stairHeight[i] = stairHeight[i - 1] + A / stairWidth[i];
    }
    return true;
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
        double x = UniformRand::StandardVariate(randGenerator) * stairWidth[stairId];
        if (x < stairWidth[stairId + 1]) /// if we are under the upper stair - accept
            return x;
        if (stairId == 0) /// if we catch the tail
            return x1 + StandardVariate(randGenerator);
        long double height = stairHeight[stairId] - stairHeight[stairId - 1];
        if (stairHeight[stairId - 1] + height * UniformRand::StandardVariate(randGenerator) < std::exp(-x)) /// if we are under the curve - accept
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
