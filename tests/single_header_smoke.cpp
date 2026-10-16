#include "randlib_single.hpp"

#include <random>
#include <vector>

int main()
{
    std::mt19937_64 engine(42);
    const randlib::NormalRand<double> normal(0.0, 1.0);
    std::vector<double> sample(32);
    normal.Sample(engine, sample);
    const auto fit = randlib::fit<randlib::NormalRand<double>>(sample);
    return fit.converged ? 0 : 1;
}
