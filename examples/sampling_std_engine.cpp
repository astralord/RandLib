#include <randlib/continuous/normal.hpp>

#include <iostream>
#include <random>
#include <vector>

int main()
{
    std::mt19937 engine(42);
    const randlib::NormalRand<double> normal(2.0, 0.5);

    std::vector<double> sample(8);
    normal.Sample(engine, sample);

    std::cout << "first draw: " << sample.front() << '\n';
}
