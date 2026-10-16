// Internal coherence of every distribution on its golden x grid:
// f = exp(logf), F + S = 1, F is monotone, Quantile inverts F and Quantile1m inverts S,
// and the density integrates to the increments of F.

#include <catch2/catch_test_macros.hpp>

#include "data.h"
#include "registry.h"
#include "stats.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

using namespace randlib_test;

namespace
{

std::map<std::string, std::vector<double>> gridByCase()
{
    std::map<std::string, std::vector<double>> grid;
    std::set<std::string> goldenCases;
    for (const Record &r : loadRecords("golden_v1")) {
        if (r.function == "F" && r.status == "ok") {
            grid[r.caseId].push_back(parseNumber(r.arg));
            goldenCases.insert(r.caseId);
        }
    }
    for (const Record &r : loadRecords("reference")) {
        if (r.function == "F" && !goldenCases.contains(r.caseId))
            grid[r.caseId].push_back(parseNumber(r.arg));
    }
    for (auto &entry : grid) {
        auto &xs = entry.second;
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    }
    return grid;
}

std::string at(const std::string &what, double x)
{
    std::ostringstream s;
    s.precision(17);
    s << what << " at x = " << x;
    return s.str();
}

std::string at(const std::string &what, double x, double actual, double expected)
{
    std::ostringstream s;
    s.precision(17);
    s << what << " at x = " << x << ": " << actual << " vs " << expected;
    return s.str();
}

void report(const std::vector<std::string> &ids, const std::vector<std::vector<std::string>> &failures)
{
    size_t total = 0;
    for (size_t i = 0; i < ids.size(); ++i) {
        for (const auto &f : failures[i]) {
            if (f.rfind("SKIPPED: ", 0) == 0) {
                WARN(ids[i] + ": " + f);
                continue;
            }
            FAIL_CHECK(ids[i] + ": " + f);
            ++total;
        }
    }
    CHECK(total == 0);
}

} // namespace

TEST_CASE("Continuous distributions are internally coherent", "[coherence]")
{
    auto grid = gridByCase();
    auto cases = casesOfKind(Kind::Continuous);
    std::vector<std::string> ids;
    for (const Case *c : cases)
        ids.push_back(c->id);
    std::vector<std::vector<std::string>> failures(cases.size());

    parallelFor(cases.size(), [&](size_t i) {
        auto &out = failures[i];
        std::unique_ptr<ContinuousDistribution<double>> d;
        try {
            d = cases[i]->makeContinuous();
        } catch (const std::exception &e) {
            out.push_back(std::string("construction failed: ") + e.what());
            return;
        }
        std::vector<double> xs = grid[cases[i]->id];
        if (xs.empty())
            return;
        bool expensive = false;
        try {
            auto start = std::chrono::steady_clock::now();
            static_cast<void>(d->F(xs.front()));
            static_cast<void>(d->f(xs.front()));
            expensive = std::chrono::steady_clock::now() - start > std::chrono::milliseconds(1);
        } catch (const randlib::not_implemented &e) {
            out.push_back(std::string("SKIPPED: ") + e.what());
            return;
        }
        double previousF = 0.0, previousX = xs.front();
        for (size_t j = 0; j < xs.size(); ++j) {
            double x = xs[j];
            /// For distributions whose f and F are themselves numerical integrals, check the
            /// quantiles at every third point and skip the quadrature of f
            bool full = !expensive || j % 3 == 0;
            try {
                double f = d->f(x), logf = d->logf(x), F = d->F(x), S = d->S(x);
                if (f > 1e-300 && std::isfinite(f) && !relativelyClose(f, std::exp(logf), 1e-12))
                    out.push_back(at("f != exp(logf)", x, f, std::exp(logf)));
                if (f == 0.0 && std::exp(logf) != 0.0)
                    out.push_back(at("f underflows before exp(logf)", x));
                if (!(std::fabs(F + S - 1.0) <= 1e-12))
                    out.push_back(at("F + S != 1", x, F + S, 1.0));
                if (F < previousF - 1e-15)
                    out.push_back(at("F decreases", x, F, previousF));
                if (full && F > 1e-6 && F < 1 - 1e-6) {
                    double q = d->Quantile(F);
                    if (!relativelyClose(q, x, 1e-7) && std::fabs(q - x) > 1e-9)
                        out.push_back(at("Quantile(F(x)) != x", x, q, x));
                }
                if (full && S > 1e-6 && S < 1 - 1e-6) {
                    double q = d->Quantile1m(S);
                    if (!relativelyClose(q, x, 1e-7) && std::fabs(q - x) > 1e-9)
                        out.push_back(at("Quantile1m(S(x)) != x", x, q, x));
                }
                /// Skip intervals that touch the support boundary (where a density may be singular
                /// beyond what doubles can resolve, or an atom may sit) and point masses
                bool touchesBoundary = previousX <= d->MinValue() || x >= d->MaxValue();
                bool pointMass = d->Variance() == 0;
                if (!expensive && !touchesBoundary && !pointMass && j > 0 && std::isfinite(previousX) &&
                    F - previousF > 1e-10) {
                    /// An integrable singularity at a support boundary has measure zero
                    double integral = integrate(
                        [&d](double t) {
                            double v = d->f(t);
                            return std::isfinite(v) ? v : 0.0;
                        },
                        previousX, x, 1e-10);
                    if (!(std::fabs(integral - (F - previousF)) <= 1e-8 + 1e-8 * (F - previousF)))
                        out.push_back(at("integral of f != F(b) - F(a)", x, integral, F - previousF));
                }
                previousF = F;
                previousX = x;
            } catch (const std::exception &e) {
                out.push_back(at(std::string("exception: ") + e.what(), x));
            }
        }
    });
    report(ids, failures);
}

TEST_CASE("Discrete distributions are internally coherent", "[coherence]")
{
    std::map<std::string, std::vector<double>> grid;
    for (const Record &r : loadRecords("golden_v1")) {
        if (r.function == "P" && r.status == "ok")
            grid[r.caseId].push_back(parseNumber(r.arg));
    }
    auto cases = casesOfKind(Kind::Discrete);
    std::vector<std::string> ids;
    for (const Case *c : cases)
        ids.push_back(c->id);
    std::vector<std::vector<std::string>> failures(cases.size());

    parallelFor(cases.size(), [&](size_t i) {
        auto &out = failures[i];
        std::unique_ptr<DiscreteDistribution<int>> d;
        try {
            d = cases[i]->makeDiscrete();
        } catch (const std::exception &e) {
            out.push_back(std::string("construction failed: ") + e.what());
            return;
        }
        auto ks = grid[cases[i]->id];
        std::sort(ks.begin(), ks.end());
        double previousF = 0.0;
        for (double kd : ks) {
            int k = static_cast<int>(kd);
            try {
                double P = d->P(k), logP = d->logP(k), F = d->F(k), S = d->S(k), Fprev = d->F(k - 1);
                if (P > 1e-300 && !relativelyClose(P, std::exp(logP), 1e-12))
                    out.push_back(at("P != exp(logP)", k, P, std::exp(logP)));
                if (P == 0.0 && logP != -INFINITY)
                    out.push_back(at("P == 0 but logP is finite", k));
                if (!(std::fabs(F + S - 1.0) <= 1e-12))
                    out.push_back(at("F + S != 1", k, F + S, 1.0));
                if (!(std::fabs(F - Fprev - P) <= 1e-12 + 1e-10 * P))
                    out.push_back(at("F(k) - F(k-1) != P(k)", k, F - Fprev, P));
                if (F < previousF - 1e-15)
                    out.push_back(at("F decreases", k, F, previousF));
                if (P > 1e-12 && F > 1e-6 && F < 1 - 1e-6) {
                    int q = d->Quantile(F);
                    if (q != k)
                        out.push_back(at("Quantile(F(k)) != k", k, q, k));
                }
                if (P > 1e-12 && Fprev > 1e-6 && S > 1e-6 && S < 1 - 1e-6) {
                    int q = d->Quantile1m(S);
                    if (q != k)
                        out.push_back(at("Quantile1m(S(k)) != k", k, q, k));
                }
                previousF = F;
            } catch (const std::exception &e) {
                out.push_back(at(std::string("exception: ") + e.what(), k));
            }
        }
    });
    report(ids, failures);
}

TEST_CASE("Singular distributions are internally coherent", "[coherence]")
{
    auto grid = gridByCase();
    for (const Case *c : casesOfKind(Kind::Singular)) {
        INFO(c->id);
        auto d = c->makeSingular();
        double previousF = 0.0;
        for (double x : grid[c->id]) {
            double F = d->F(x), S = d->S(x);
            CHECK(std::fabs(F + S - 1.0) <= 1e-12);
            CHECK(F >= previousF);
            previousF = F;
        }
    }
}
