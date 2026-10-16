#include <randlib/randlib.hpp>

#include <iostream>
#include <vector>

int main()
{
    const std::vector<double> observations = {1.2, 0.7, 1.8, 1.1, 0.9};
    const randlib::NormalInverseGammaRand<double> prior(0.0, 1.0, 2.0, 2.0);
    const auto posterior = randlib::fit_bayes(observations, prior);

    std::cout << "posterior mean location: " << posterior.GetLocation() << '\n'
              << "posterior precision: " << posterior.GetPrecision() << '\n';
}
