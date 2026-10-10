// Golden regression: the library must reproduce the values recorded from the v1.1
// baseline (tests/data/golden_v1, written by tools/golden/run_dump.sh), except where
// tests/data/golden_v1/known_deviations.csv records an intentional change.

#include <catch2/catch_test_macros.hpp>

#include "data.h"
#include "registry.h"
#include "stats.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

using namespace randlib_test;

namespace
{

struct Deviation
{
    std::string caseId, function, arg, expected, reason;

    bool matches(const Record &r) const
    {
        return (caseId == "*" || caseId == r.caseId) && (function == "*" || function == r.function)
               && (arg == "*" || arg == r.arg);
    }
};

std::vector<Deviation> loadDeviations()
{
    std::ifstream in(dataDirectory() + "/golden_v1/known_deviations.csv");
    std::vector<Deviation> deviations;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        auto f = splitCsvLine(line);
        deviations.push_back({f.at(0), f.at(1), f.at(2), f.at(3), f.at(4)});
    }
    return deviations;
}

struct Tolerance
{
    double relative;
    double absolute;
};

Tolerance goldenTolerance(const std::string &function)
{
    if (function == "MinValue" || function == "MaxValue")
        return {0.0, 0.0};
    if (function == "Quantile" || function == "Quantile1m")
        return {1e-7, 1e-12};
    if (function == "Mode" || function == "Median" || function == "Mode1" || function == "Mode2")
        return {1e-6, 1e-9};
    if (function == "logf" || function == "logP" || function.rfind("CF", 0) == 0)
        return {1e-9, 1e-12};
    return {1e-9, 0.0};
}

bool close(double actual, double expected, Tolerance t)
{
    if (relativelyClose(actual, expected, t.relative))
        return true;
    return std::isfinite(actual) && std::isfinite(expected) && std::fabs(actual - expected) <= t.absolute;
}

std::string describe(const Record &r)
{
    return r.caseId + " " + r.function + "(" + r.arg + ")";
}

/// Returns an empty string when the record is reproduced, otherwise a failure message.
/// The actual outcome is stored in `actualText` for the optional mismatch report.
std::string check(const Record &r, const Evaluator &evaluate, const std::vector<Deviation> &deviations, std::string &actualText)
{
    std::string expected = r.status == "ok" ? r.value : r.status;
    for (const Deviation &d : deviations) {
        if (d.matches(r)) {
            if (d.expected == "skip")
                return {};
            expected = d.expected;
            break;
        }
    }

    std::string outcome;
    double actual = 0.0;
    try {
        actual = evaluate(r.function, r.arg);
        outcome = "ok";
        std::ostringstream text;
        text.precision(17);
        text << actual;
        actualText = text.str();
    } catch (const std::invalid_argument &) {
        outcome = "throw:invalid_argument";
    } catch (const std::runtime_error &) {
        outcome = "throw:runtime_error";
    } catch (const std::exception &e) {
        outcome = std::string("throw:exception (") + e.what() + ")";
    }

    if (outcome != "ok")
        actualText = outcome;
    if (expected.rfind("throw", 0) == 0)
        return outcome == expected ? std::string{} : describe(r) + ": expected " + expected + ", got " + outcome;
    if (outcome != "ok")
        return describe(r) + ": expected " + expected + ", got " + outcome;
    if (expected == "no-throw")
        return {};

    Record wanted = r;
    wanted.value = expected;
    if (close(actual, wanted.number(), goldenTolerance(r.function)))
        return {};
    std::ostringstream message;
    message.precision(17);
    message << describe(r) << ": expected " << expected << ", got " << actual
            << " (relative error " << relativeError(actual, wanted.number()) << ")";
    return message.str();
}

}

TEST_CASE("Golden v1 values are reproduced", "[golden]")
{
    const auto records = loadRecords("golden_v1");
    const auto deviations = loadDeviations();
    REQUIRE(records.size() > 10000);

    std::map<std::string, std::vector<const Record *>> byCase;
    for (const Record &r : records) {
        if (r.function == "Name" || r.status == "timeout")
            continue;
        byCase[r.caseId].push_back(&r);
    }
    std::vector<std::string> caseIds;
    for (const auto &entry : byCase)
        caseIds.push_back(entry.first);

    std::vector<std::vector<std::string>> failures(caseIds.size());
    std::vector<std::vector<std::string>> reportRows(caseIds.size());
    parallelFor(caseIds.size(), [&](size_t i) {
        const Case *c = findCase(caseIds[i]);
        if (c == nullptr) {
            failures[i].push_back(caseIds[i] + ": no such case in tests/support/cases.inc");
            return;
        }
        Evaluator evaluate;
        std::string constructed = "ok";
        try {
            evaluate = c->makeEvaluator();
        } catch (const std::invalid_argument &) {
            constructed = "throw:invalid_argument";
        }
        for (const Record *r : byCase[caseIds[i]]) {
            if (r->function == "__ctor__") {
                std::string expected = r->status;
                for (const Deviation &d : deviations) {
                    if (d.matches(*r))
                        expected = d.expected == "no-throw" ? "ok" : d.expected;
                }
                if (constructed != expected)
                    failures[i].push_back(describe(*r) + ": expected " + expected + ", got " + constructed);
                continue;
            }
            if (!evaluate)
                continue;
            std::string actualText;
            std::string failure = check(*r, evaluate, deviations, actualText);
            if (!failure.empty()) {
                failures[i].push_back(failure);
                reportRows[i].push_back("\"" + r->caseId + "\"," + r->function + "," + r->arg + "," + r->value + "," + r->status + "," + actualText);
            }
        }
    });

    /// RANDLIB_GOLDEN_REPORT=<file> writes the mismatches as CSV, for tools/golden/classify_deviations.py
    if (const char *reportPath = std::getenv("RANDLIB_GOLDEN_REPORT")) {
        std::ofstream report(reportPath);
        report << "case,function,arg,golden,status,actual\n";
        for (const auto &rows : reportRows) {
            for (const auto &row : rows)
                report << row << "\n";
        }
    }

    size_t total = 0;
    for (const auto &caseFailures : failures) {
        for (const auto &failure : caseFailures) {
            ++total;
            FAIL_CHECK(failure);
        }
    }
    CHECK(total == 0);
}
