#include <randlib/randlib.hpp>

#include <iostream>
#include <string>
#include <variant>
#include <vector>

int main()
{
    std::vector<double> sample;
    sample.reserve(400);
    const randlib::ExponentialRand<double> reference(2.0);
    for (int i = 0; i < 400; ++i)
        sample.push_back(reference.Quantile((i + 0.5) / 400));

    const auto best =
        randlib::fit_best<randlib::ExponentialRand<double>, randlib::GammaRand<double>, randlib::LogNormalRand<double>>(
            sample, randlib::Criterion::AIC);
    const std::string name = std::visit([](const auto &fit) { return fit.distribution.Name(); }, best.fit);
    std::cout << "best model: " << name << " (AIC " << best.value << ")\n";
}
