#include <randlib/randlib.hpp>

#include <iostream>
#include <vector>

int main()
{
    using Stable = randlib::StableRand<double>;
    const Stable source(1.4, 0.35, 1.5, -0.5, Stable::Parameterization::S0);
    randlib::Pcg64 engine(42);
    std::vector<double> sample(2000);
    source.Sample(engine, sample);

    const auto estimate = randlib::fit<Stable>(sample, randlib::Method::Quantile);
    const auto &model = estimate.distribution;
    std::cout << "alpha=" << model.GetExponent() << " beta=" << model.GetSkewness() << " scale=" << model.GetScale()
              << " location=" << model.GetLocation() << '\n';
}
