// Records the numerical behaviour of the v1.1 baseline library as CSV.
//
//   dump_v1 --list        print "<index>\t<case id>" for every case
//   dump_v1 --case N [--skip KEY,KEY...]
//                         print CSV rows for case N to stdout; each KEY is
//                         "function@arg" and marks a call that hangs in the baseline
//
// Build and run through tools/golden/run_dump.sh, which links this file against the
// library at tag v1.1.0-baseline and runs every case in its own process with a timeout.

#include "RandLib.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace {

std::string caseId(std::string type, const std::string &args)
{
    const std::string plain = "RL_PLAIN(";
    if (type.rfind(plain, 0) == 0)
        type = type.substr(plain.size(), type.size() - plain.size() - 1);
    auto lt = type.find('<');
    if (lt != std::string::npos)
        type.resize(lt);
    return type + "(" + args + ")";
}

std::string fmt(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

std::string csvQuote(const std::string &s)
{
    std::string out = "\"";
    for (char c : s) {
        if (c == '"')
            out += '"';
        out += c;
    }
    return out + "\"";
}

struct Out
{
    std::string id;
    std::vector<std::string> skipped;

    void row(const std::string &function, const std::string &arg, const std::string &value, const char *status) const
    {
        std::printf("%s,%s,%s,%s,%s\n", csvQuote(id).c_str(), function.c_str(), arg.c_str(), value.c_str(), status);
        std::fflush(stdout);
    }

    void value(const std::string &function, const std::string &arg, const std::function<double()> &fn, bool record = true) const
    {
        const std::string key = function + "@" + arg;
        if (std::find(skipped.begin(), skipped.end(), key) != skipped.end()) {
            row(function, arg, "", "timeout");
            return;
        }
        std::fprintf(stderr, "CALL %s\n", key.c_str());
        std::fflush(stderr);
        try {
            double result = fn();
            if (record)
                row(function, arg, fmt(result), "ok");
        } catch (const std::invalid_argument &) {
            row(function, arg, "", "throw:invalid_argument");
        } catch (const std::runtime_error &) {
            row(function, arg, "", "throw:runtime_error");
        } catch (const std::exception &) {
            row(function, arg, "", "throw:exception");
        }
    }

    void complexValue(const std::string &function, const std::string &arg, const std::function<std::complex<double>()> &fn) const
    {
        std::complex<double> z;
        bool ok = false;
        value(function, arg, [&] { z = fn(); ok = true; return 0.0; }, false);
        if (!ok)
            return;
        row(function + "_re", arg, fmt(z.real()), "ok");
        row(function + "_im", arg, fmt(z.imag()), "ok");
    }
};

double roundSignificant(double x, int digits)
{
    if (x == 0.0 || !std::isfinite(x))
        return x;
    double scale = std::pow(10.0, digits - 1 - static_cast<int>(std::floor(std::log10(std::fabs(x)))));
    return std::round(x * scale) / scale;
}

const std::vector<double> SAMPLE_LEVELS = {0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999};
const std::vector<double> P_GRID = {1e-10, 1e-6, 0.001, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999, 1 - 1e-6};
const std::vector<double> T_GRID = {0.1, 0.5, 1, 2, 5};
constexpr unsigned long GRID_SEED = 12345;
constexpr size_t GRID_SAMPLE_SIZE = 20000;

template <class D>
std::vector<double> continuousGrid(const D &d)
{
    std::vector<double> grid;
    std::vector<double> sample(GRID_SAMPLE_SIZE);
    try {
        d.Reseed(GRID_SEED);
        d.Sample(sample);
        std::sort(sample.begin(), sample.end());
        for (double level : SAMPLE_LEVELS)
            grid.push_back(roundSignificant(sample[static_cast<size_t>(level * (sample.size() - 1))], 4));
    } catch (const std::exception &) {
    }
    double minV = d.MinValue(), maxV = d.MaxValue();
    if (std::isfinite(minV)) {
        grid.push_back(minV);
        grid.push_back(minV - 1);
    }
    if (std::isfinite(maxV)) {
        grid.push_back(maxV);
        grid.push_back(maxV + 1);
    }
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
    return grid;
}

template <class D>
std::vector<int> discreteGrid(const D &d)
{
    std::vector<int> grid;
    std::vector<int> sample(GRID_SAMPLE_SIZE);
    try {
        d.Reseed(GRID_SEED);
        d.Sample(sample);
        std::sort(sample.begin(), sample.end());
        for (double level : SAMPLE_LEVELS) {
            int k = sample[static_cast<size_t>(level * (sample.size() - 1))];
            grid.push_back(k - 1);
            grid.push_back(k);
            grid.push_back(k + 1);
        }
    } catch (const std::exception &) {
    }
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
    return grid;
}

template <class D>
void dumpMoments(const Out &o, const D &d)
{
    o.row("Name", "", csvQuote(d.Name()), "ok");
    o.value("MinValue", "", [&] { return static_cast<double>(d.MinValue()); });
    o.value("MaxValue", "", [&] { return static_cast<double>(d.MaxValue()); });
    o.value("Mean", "", [&] { return static_cast<double>(d.Mean()); });
    o.value("Variance", "", [&] { return static_cast<double>(d.Variance()); });
    o.value("Skewness", "", [&] { return static_cast<double>(d.Skewness()); });
    o.value("ExcessKurtosis", "", [&] { return static_cast<double>(d.ExcessKurtosis()); });
    o.value("Median", "", [&] { return static_cast<double>(d.Median()); });
}

template <class D>
void dumpQuantilesAndCF(const Out &o, const D &d)
{
    for (double p : P_GRID) {
        o.value("Quantile", fmt(p), [&] { return static_cast<double>(d.Quantile(p)); });
        o.value("Quantile1m", fmt(p), [&] { return static_cast<double>(d.Quantile1m(p)); });
    }
    for (double t : T_GRID)
        o.complexValue("CF", fmt(t), [&] { return d.CF(t); });
}

void dumpContinuous(const Out &o, const ContinuousDistribution<double> &d)
{
    dumpMoments(o, d);
    o.value("Mode", "", [&] { return static_cast<double>(d.Mode()); });
    for (double x : continuousGrid(d)) {
        o.value("f", fmt(x), [&] { return d.f(x); });
        o.value("logf", fmt(x), [&] { return d.logf(x); });
        o.value("F", fmt(x), [&] { return d.F(x); });
        o.value("S", fmt(x), [&] { return d.S(x); });
    }
    dumpQuantilesAndCF(o, d);
}

void dumpSingular(const Out &o, const UnivariateDistribution<double> &d)
{
    dumpMoments(o, d);
    for (double x : continuousGrid(d)) {
        o.value("F", fmt(x), [&] { return d.F(x); });
        o.value("S", fmt(x), [&] { return d.S(x); });
    }
    dumpQuantilesAndCF(o, d);
}

void dumpDiscrete(const Out &o, const DiscreteDistribution<int> &d)
{
    dumpMoments(o, d);
    o.value("Mode", "", [&] { return static_cast<double>(d.Mode()); });
    for (int k : discreteGrid(d)) {
        std::string arg = std::to_string(k);
        o.value("P", arg, [&] { return d.P(k); });
        o.value("logP", arg, [&] { return d.logP(k); });
        o.value("F", arg, [&] { return d.F(k); });
        o.value("S", arg, [&] { return d.S(k); });
    }
    dumpQuantilesAndCF(o, d);
}

template <class D>
void dumpBivariateMoments(const Out &o, const D &d)
{
    o.row("Name", "", csvQuote(d.Name()), "ok");
    o.value("Mean1", "", [&] { return static_cast<double>(d.Mean().first); });
    o.value("Mean2", "", [&] { return static_cast<double>(d.Mean().second); });
    o.value("Cov11", "", [&] { return static_cast<double>(std::get<0>(d.Covariance())); });
    o.value("Cov12", "", [&] { return static_cast<double>(std::get<1>(d.Covariance())); });
    o.value("Cov22", "", [&] { return static_cast<double>(std::get<2>(d.Covariance())); });
    o.value("Correlation", "", [&] { return static_cast<double>(d.Correlation()); });
    o.value("Mode1", "", [&] { return static_cast<double>(d.Mode().first); });
    o.value("Mode2", "", [&] { return static_cast<double>(d.Mode().second); });
}

template <class D>
void dumpBivariateContinuous(const Out &o, const D &d)
{
    dumpBivariateMoments(o, d);
    const std::vector<double> axis = {-2.0, -0.5, 0.25, 1.0, 3.0};
    for (double x : axis) {
        for (double y : axis) {
            std::string arg = fmt(x) + ";" + fmt(y);
            o.value("f", arg, [&] { return d.f(DoublePair(x, y)); });
            o.value("logf", arg, [&] { return d.logf(DoublePair(x, y)); });
        }
    }
}

template <class D>
void dumpBivariateDiscrete(const Out &o, const D &d)
{
    dumpBivariateMoments(o, d);
    for (int x = -1; x <= 6; ++x) {
        for (int y = -1; y <= 6; ++y) {
            std::string arg = std::to_string(x) + ";" + std::to_string(y);
            o.value("P", arg, [&] { return d.P(IntPair(x, y)); });
            o.value("logP", arg, [&] { return d.logP(IntPair(x, y)); });
        }
    }
}

struct Case
{
    std::string id;
    std::function<void(const Out &)> run;
};

#define RL_PLAIN(X) X
#define RL_ADD(dumper, T, ...) \
    cases.push_back({caseId(#T, #__VA_ARGS__), [](const Out &o) { \
        std::unique_ptr<const T> d; \
        o.value("__ctor__", "", [&] { d = std::make_unique<const T>(__VA_ARGS__); return 0.0; }); \
        if (d) \
            dumper(o, *d); \
    }});
#define RL_CONT(T, ...) RL_ADD(dumpContinuous, T, __VA_ARGS__)
#define RL_CONT_V2(T, ...)
#define RL_CONT_NO_STATS(T, ...) RL_CONT(T, __VA_ARGS__)
#define RL_CONT_V2_NO_STATS(T, ...)
#define RL_DISC(T, ...) RL_ADD(dumpDiscrete, T, __VA_ARGS__)
#define RL_SING(T, ...) RL_ADD(dumpSingular, T, __VA_ARGS__)
#define RL_BIVC(T, ...) RL_ADD(dumpBivariateContinuous, T, __VA_ARGS__)
#define RL_BIVD(T, ...) RL_ADD(dumpBivariateDiscrete, T, __VA_ARGS__)

std::vector<Case> allCases()
{
    std::vector<Case> cases;
#include "../../tests/support/cases.inc"
    return cases;
}

} // namespace

int main(int argc, char **argv)
{
    const auto cases = allCases();
    if (argc == 2 && std::strcmp(argv[1], "--list") == 0) {
        for (size_t i = 0; i != cases.size(); ++i)
            std::printf("%zu\t%s\n", i, cases[i].id.c_str());
        return 0;
    }
    if ((argc == 3 || argc == 5) && std::strcmp(argv[1], "--case") == 0) {
        size_t index = std::stoul(argv[2]);
        if (index >= cases.size())
            return 2;
        Out out{cases[index].id, {}};
        if (argc == 5 && std::strcmp(argv[3], "--skip") == 0) {
            std::string list = argv[4];
            size_t start = 0;
            while (start <= list.size()) {
                size_t comma = list.find(',', start);
                if (comma == std::string::npos)
                    comma = list.size();
                if (comma > start)
                    out.skipped.push_back(list.substr(start, comma - start));
                start = comma + 1;
            }
        }
        cases[index].run(out);
        return 0;
    }
    std::fprintf(stderr, "usage: %s --list | --case N\n", argv[0]);
    return 2;
}
