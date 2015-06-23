#include "ExponentialRand.h"

ExponentialRand::ExponentialRand(double rate) :
    U(0, 1)
{
    setRate(rate);

    /// Set up ziggurat tables
    const double m2 = 4294967296.;
    double de = 7.697117470131487, te = de, ve = 3.949659822581572e-3;

    fe[0] = 1.;
    fe[255] = std::exp(-de);

    double q = ve / fe[255];
    ke[0] = (de / q) * m2;
    ke[1] = 0;

    we[0] = q / m2;
    we[255] = de / m2;

    for (size_t i = 254; i >= 1; --i)
    {
       de = -std::log(ve / de + std::exp(-de));
       ke[i + 1] = (de / te) * m2;
       te = de;
       fe[i] = std::exp(-de);
       we[i] = de / m2;
    }
}

void ExponentialRand::setRate(double rate)
{
    l = std::max(rate, MIN_POSITIVE);
    beta = 1.0 / l;
}

double ExponentialRand::f(double x) const
{
    return (x > 0) ? l * std::exp(-l * x) : 0;
}

double ExponentialRand::F(double x) const
{
    return (x > 0) ? 1 - std::exp(-l * x) : 0;
}

double ExponentialRand::ziggurat()
{
    for (;;)
    {
        unsigned long jz = SHR3();
        unsigned long iz = jz & 255;
        double x = jz * we[iz];

        if (jz < ke[iz])
            return x;

        if (iz == 0)
            return (7.69711 - std::log(U.value()));

        if (fe[iz] + U.value() * (fe[iz - 1] - fe[iz]) < std::exp(-x))
            return x;
    }
}

double ExponentialRand::value()
{
    return beta * ziggurat();
}
