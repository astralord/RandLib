// Mode() is a local maximum of the density (or mass function) and Median() splits the mass in half.

#include <catch2/catch_test_macros.hpp>

#include "registry.h"
#include "stats.h"

#include <cmath>
#include <sstream>

using namespace randlib_test;

TEST_CASE("Continuous mode maximizes the density and median halves the mass", "[mode]")
{
    auto cases = casesOfKind(Kind::Continuous);
    std::vector<std::string> failures(cases.size());
    parallelFor(cases.size(), [&](size_t i) {
        std::ostringstream out;
        out.precision(17);
        try {
            auto d = cases[i]->makeContinuous();
            double median = d->Median();
            double Fmedian = d->F(median);
            bool pointMass = d->Variance() == 0;
            if (!pointMass && !(std::fabs(Fmedian - 0.5) <= 1e-8))
                out << "F(Median()) = " << Fmedian << " at median " << median << "; ";

            double mode = d->Mode();
            double scale = d->Quantile(0.75) - d->Quantile(0.25);
            double h = 1e-3 * scale;
            double fMode = d->f(mode);
            for (double x : {mode - h, mode + h}) {
                if (x < d->MinValue() || x > d->MaxValue())
                    continue;
                double fx = d->f(x);
                if (fx > fMode * (1 + 1e-9) && std::isfinite(fMode))
                    out << "f(" << x << ") = " << fx << " > f(Mode() = " << mode << ") = " << fMode << "; ";
            }
        } catch (const randlib::not_implemented &e) {
            out << "SKIPPED: " << e.what();
        } catch (const std::exception &e) {
            out << "exception: " << e.what();
        }
        failures[i] = out.str();
    });
    for (size_t i = 0; i < cases.size(); ++i) {
        if (failures[i].rfind("SKIPPED: ", 0) == 0)
            WARN(cases[i]->id + ": " + failures[i]);
        else if (!failures[i].empty())
            FAIL_CHECK(cases[i]->id + ": " + failures[i]);
    }
}

TEST_CASE("Discrete mode maximizes the mass function and median halves the mass", "[mode]")
{
    for (const Case *c : casesOfKind(Kind::Discrete)) {
        INFO(c->id);
        std::unique_ptr<DiscreteDistribution<int>> d;
        try {
            d = c->makeDiscrete();
        } catch (const std::exception &e) {
            FAIL_CHECK(c->id + ": construction failed: " + e.what());
            continue;
        }
        int median = d->Median();
        CHECK(d->F(median) >= 0.5 - 1e-12);
        CHECK(d->F(median - 1) <= 0.5 + 1e-12);
        int mode = d->Mode();
        double pMode = d->P(mode);
        CHECK(d->P(mode - 1) <= pMode * (1 + 1e-12));
        CHECK(d->P(mode + 1) <= pMode * (1 + 1e-12));
    }
}
