#include <randlib/randlib.hpp>

#include <cstddef>
#include <iostream>
#include <thread>
#include <vector>

int main()
{
    auto streams = randlib::make_streams(42, 4);
    const randlib::ExponentialRand<double> exponential(2.0);
    std::vector<double> means(streams.size());
    std::vector<std::jthread> workers;
    workers.reserve(streams.size());

    for (std::size_t i = 0; i < streams.size(); ++i) {
        workers.emplace_back([&, i] {
            double sum = 0;
            for (int draw = 0; draw < 10'000; ++draw)
                sum += exponential(streams[i]);
            means[i] = sum / 10'000;
        });
    }
    workers.clear();

    for (std::size_t i = 0; i < means.size(); ++i)
        std::cout << "stream " << i << " mean: " << means[i] << '\n';
}
