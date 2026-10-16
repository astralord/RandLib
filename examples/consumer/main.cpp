// Builds against an installed RandLib through find_package(RandLib); used by CI.
#include <randlib/randlib.hpp>

#include <cstdio>
#include <random>
#include <vector>

int main()
{
    randlib::NormalRand<double> X(0, 1);
    std::mt19937_64 engine(1);
    std::vector<double> sample(1000);
    X.Sample(engine, sample);
    std::vector<double> gammaSample(1000);
    randlib::GammaRand<double>(2.5, 0.5).Sample(engine, gammaSample);
    auto fitted = randlib::fit<randlib::GammaRand<double>>(gammaSample).distribution;
    std::printf("%s, F(0) = %.6f, fitted %s\n", X.Name().c_str(), X.F(0.0), fitted.Name().c_str());
    return (X.F(0.0) == 0.5) ? 0 : 1;
}
