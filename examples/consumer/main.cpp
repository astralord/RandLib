// Builds against an installed RandLib through find_package(RandLib); used by CI.
#include <RandLib.h>

#include <cstdio>
#include <vector>

int main()
{
    NormalRand<double> X(0, 1);
    X.Reseed(1);
    std::vector<double> sample(1000);
    X.Sample(sample);
    GammaRand<double> fitted;
    std::vector<double> gammaSample(1000);
    GammaRand<double>(2.5, 0.5).Sample(gammaSample);
    fitted.Fit(gammaSample);
    std::printf("%s, F(0) = %.6f, fitted %s\n", X.Name().c_str(), X.F(0.0), fitted.Name().c_str());
    return (X.F(0.0) == 0.5) ? 0 : 1;
}
