# Parallel sampling and independent streams

Distribution objects contain no shared mutable sampling state. Every `operator()(engine)` and
`Sample(engine, output)` call advances only the engine supplied by the caller, so parallel ownership
is explicit.

```cpp
#include <randlib/randlib.hpp>

#include <thread>
#include <vector>

auto streams = randlib::make_streams(42, 8);
const randlib::NormalRand<double> normal(0.0, 1.0);
std::vector<std::vector<double>> samples(8, std::vector<double>(100'000));
std::vector<std::jthread> workers;

for (std::size_t i = 0; i < streams.size(); ++i) {
    workers.emplace_back([&, i] { normal.Sample(streams[i], samples[i]); });
}
```

Keep one engine per worker and do not use the same engine concurrently. `make_streams` chooses the
separation strategy from the engine type:

- `Xoshiro256pp` (the default) uses successive `jump()` operations, separating streams by 2^128
  states;
- `Pcg32` and `Pcg64` use distinct PCG sequence parameters; and
- other standard engines, including `std::mt19937_64`, receive seeds expanded by `SplitMix64`.

The first two strategies have structural stream separation. Derived seeds for arbitrary engines are
convenient and deterministic, but they do not prove non-overlap.

For reproducible work, record the engine type, seed, stream index, compiler/platform, and RandLib
version. The C++ standard specifies distribution results statistically rather than requiring the same
sample bits across implementations. RandLib likewise reserves the right to fix sampling defects;
v1.2 and v2 intentionally do not reproduce v1.1 streams.

Use `random_seed()` when reproducibility is not wanted. It mixes `std::random_device`, a high-resolution
clock, the thread id, and `SplitMix64`; it is not a cryptographic random-number facility.
