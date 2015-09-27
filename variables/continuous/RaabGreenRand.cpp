#include "RaabGreenRand.h"
#include "UniformRand.h"

RaabGreenRand::RaabGreenRand()
{

}

std::string RaabGreenRand::name()
{
    return "Raab-Green";
}

double RaabGreenRand::f(double x) const
{
    double absX = std::fabs(x);
    return absX < M_PI ? (0.5 * (1.0 + std::cos(x)) * M_1_PI) : 0.0;
}

double RaabGreenRand::F(double x) const
{
    if (x <= -M_PI)
        return 0.0;
    if (x >= M_PI)
        return 1.0;
    double y = x + std::sin(x);
    y *= M_1_PI;
    ++y;
    return 0.5 * y;
}

double RaabGreenRand::variate() const
{
    /// p. 160. Non-Uniform Random Variate Generation. Luc Devroye
    double X = UniformRand::variate(-M_PI_2, M_PI_2);
    double XSq = X * X;
    double U = UniformRand::standardVariate();
    U += U;
    int a = 0, b = -1;
    double W = 0.0, V = 1.0;
    int iter = 0;
    do {
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W += V;
        if (U >= W)
            return X;

        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W -= V;
        if (U <= W)
        {
            if (X == 0.0)
                return 0.0;
            return (X > 0.0) ? M_PI - X : -M_PI - X;
        }
    } while (++iter < 1e9);
    return NAN;
}

double RaabGreenRand::Mean() const
{
    return 0.0;
}

double RaabGreenRand::Variance() const
{
    //TODO!!
    return -1.0;
}

double RaabGreenRand::Median() const
{
    return 0.0;
}

double RaabGreenRand::Mode() const
{
    return 0.0;
}

double RaabGreenRand::Skewness() const
{
    return 0.0;
}

