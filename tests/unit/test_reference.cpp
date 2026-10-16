// Accuracy against 40-digit mpmath reference values (tests/data/reference, written by
// tools/refdata/generate.py), with the tolerances of docs/MODERNIZATION_PLAN.md section 9.
// Loosened tolerances live in tests/data/reference_tolerances.csv, one row per gap, and
// each row must link the issue that tracks it.

#include <catch2/catch_test_macros.hpp>

#include "data.h"
#include "registry.h"
#include "stats.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

using namespace randlib_test;

namespace
{

const std::set<std::string> TIER1 = {"NormalRand",
                                     "LogNormalRand",
                                     "ExponentialRand",
                                     "GammaRand",
                                     "ChiSquaredRand",
                                     "ErlangRand",
                                     "BetaRand",
                                     "BetaPrimeRand",
                                     "StudentTRand",
                                     "NoncentralTRand",
                                     "CauchyRand",
                                     "LevyRand",
                                     "StableRand",
                                     "HoltsmarkRand",
                                     "LandauRand",
                                     "LogStableRand",
                                     "GeometricStableRand",
                                     "ParetoRand",
                                     "WeibullRand",
                                     "FrechetRand",
                                     "GumbelRand",
                                     "GeneralizedParetoRand",
                                     "GeneralizedExtremeValueRand",
                                     "LaplaceRand",
                                     "AsymmetricLaplaceRand",
                                     "LogisticRand",
                                     "InverseGammaRand",
                                     "InverseGaussianRand",
                                     "UniformRand",
                                     "FisherFRand",
                                     "NoncentralChiSquaredRand",
                                     "PoissonRand",
                                     "BinomialRand",
                                     "NegativeBinomialRand",
                                     "PascalRand",
                                     "GeometricRand",
                                     "BernoulliRand",
                                     "UniformDiscreteRand"};

const std::set<std::string> STABLE = {"StableRand", "HoltsmarkRand", "LandauRand", "LogStableRand"};

struct Override
{
    std::string casePattern, function, region;
    double tolerance;
    std::string issue;

    bool matches(const std::string &caseId, const std::string &fn, const std::string &rgn) const
    {
        bool caseMatches =
            casePattern == "*" ||
            (casePattern.back() == '*' ? caseId.rfind(casePattern.substr(0, casePattern.size() - 1), 0) == 0
                                       : caseId == casePattern);
        return caseMatches && (function == "*" || function == fn) && (region == "*" || region == rgn);
    }
};

std::vector<Override> loadOverrides()
{
    std::ifstream in(dataDirectory() + "/reference_tolerances.csv");
    std::vector<Override> overrides;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        auto f = splitCsvLine(line);
        overrides.push_back({f.at(0), f.at(1), f.at(2), parseNumber(f.at(3)), f.at(4)});
        REQUIRE(f.at(4).rfind("https://", 0) == 0);
    }
    return overrides;
}

struct Group
{
    double worst = 0.0;
    std::string where;
    double tolerance = 0.0;
    std::string kind;
};

/// Tolerance from section 9 of the plan; "relative" or "absolute" in kind.
double defaultTolerance(bool tier1, bool stable, const std::string &function, const std::string &region,
                        std::string &kind)
{
    kind = "relative";
    if (function == "Mean" || function == "Variance")
        return 1e-13;
    if (function == "Skewness" || function == "ExcessKurtosis")
        return 1e-11;
    if (function == "Quantile" || function == "Quantile1m" || function == "Median" || function == "Mode")
        return tier1 ? 1e-10 : 1e-7;
    if (function == "CF_re" || function == "CF_im") {
        kind = "absolute";
        return tier1 ? 1e-8 : 1e-6;
    }
    if (stable && (function == "f" || function == "F" || function == "S" || function == "logf")) {
        /// plan section 9: 1e-9 absolute for the stable density and distribution function
        kind = "absolute";
        return 1e-9;
    }
    if (function == "logf" || function == "logP")
        kind = "absolute";
    if (region == "bulk")
        return tier1 ? 1e-12 : 1e-8;
    return tier1 ? 1e-9 : 1e-6;
}

double error(double actual, double expected, const std::string &kind, const std::string &function)
{
    if (actual == expected || (std::isnan(actual) && std::isnan(expected)))
        return 0.0;
    if (!std::isfinite(actual) || !std::isfinite(expected))
        return INFINITY;
    if (kind == "absolute")
        return std::fabs(actual - expected);
    /// A reference value below the double range can only be represented as 0 or a subnormal
    if (std::fabs(expected) < DBL_MIN && std::fabs(actual) < DBL_MIN)
        return 0.0;
    /// Quantiles near 0 (the median of a symmetric law) are compared on an absolute scale
    if ((function == "Quantile" || function == "Quantile1m" || function == "Median" || function == "Mode") &&
        std::fabs(actual - expected) <= 1e-14)
        return 0.0;
    return relativeError(actual, expected);
}

void runReference(bool tier1Only)
{
    const auto records = loadRecords("reference");
    const auto overrides = loadOverrides();
    REQUIRE(!records.empty());

    std::map<std::string, std::vector<const Record *>> byCase;
    std::map<std::string, std::map<std::string, std::pair<double, double>>> tails;
    for (const Record &r : records) {
        byCase[r.caseId].push_back(&r);
        if (r.function == "F")
            tails[r.caseId][r.arg].first = r.number();
        if (r.function == "S")
            tails[r.caseId][r.arg].second = r.number();
    }
    std::vector<std::string> caseIds;
    for (const auto &entry : byCase) {
        std::string cls = entry.first.substr(0, entry.first.find('('));
        if ((TIER1.count(cls) > 0) == tier1Only)
            caseIds.push_back(entry.first);
    }

    std::vector<std::map<std::string, Group>> groups(caseIds.size());
    std::vector<std::string> missing(caseIds.size());
    parallelFor(caseIds.size(), [&](size_t i) {
        const std::string &id = caseIds[i];
        const Case *c = findCase(id);
        if (c == nullptr) {
            missing[i] = id;
            return;
        }
        std::string cls = c->className;
        bool tier1 = TIER1.count(cls) > 0, stable = STABLE.count(cls) > 0;
        Evaluator evaluate = c->makeEvaluator();
        for (const Record *r : byCase[id]) {
            std::string region = "-";
            auto tail = tails[id].find(r->arg);
            if (!r->arg.empty() && tail != tails[id].end() && r->function != "Quantile" &&
                r->function != "Quantile1m" && r->function.rfind("CF", 0) != 0)
                region = std::min(tail->second.first, tail->second.second) > 1e-6 ? "bulk" : "tail";
            std::string kind;
            double tolerance = defaultTolerance(tier1, stable, r->function, region, kind);
            for (const Override &o : overrides) {
                if (o.matches(id, r->function, region)) {
                    tolerance = o.tolerance;
                    break;
                }
            }
            double actual;
            try {
                actual = evaluate(r->function, r->arg);
            } catch (const std::exception &) {
                actual = NAN;
            }
            double err = error(actual, r->number(), kind, r->function);
            if (r->function == "logf" || r->function == "logP")
                err = std::max(0.0, err - 1e-15 * std::fabs(r->number()));
            Group &g = groups[i][r->function + " " + region];
            g.tolerance = tolerance;
            g.kind = kind;
            if (err > g.worst || (std::isnan(err) && !std::isnan(g.worst))) {
                std::ostringstream where;
                where.precision(17);
                where << r->function << "(" << r->arg << ") = " << actual << ", reference " << r->value;
                g.worst = err;
                g.where = where.str();
            }
        }
    });

    for (size_t i = 0; i < caseIds.size(); ++i) {
        if (!missing[i].empty())
            FAIL_CHECK(missing[i] + ": reference data for a case that tests/support/cases.inc does not define");
        for (const auto &[name, g] : groups[i]) {
            bool ok = g.worst <= g.tolerance;
            if (!ok) {
                std::ostringstream message;
                message << caseIds[i] << " " << name << ": worst " << g.kind << " error " << g.worst << " exceeds "
                        << g.tolerance << " at " << g.where;
                FAIL_CHECK(message.str());
            }
        }
    }
}

} // namespace

TEST_CASE("Tier 1 distributions match high-precision reference values", "[reference]")
{
    runReference(true);
}

TEST_CASE("Tier 2 distributions match high-precision reference values", "[reference]")
{
    runReference(false);
}
