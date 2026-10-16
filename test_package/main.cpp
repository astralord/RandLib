#include <randlib/randlib.hpp>

#include <random>

int main()
{
    std::mt19937 engine(7);
    const randlib::ExponentialRand<double> distribution(2.0);
    const double draw = distribution(engine);
    return draw > 0.0 ? 0 : 1;
}
