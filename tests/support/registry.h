#ifndef RANDLIB_TESTS_REGISTRY_H
#define RANDLIB_TESTS_REGISTRY_H

#include "RandLib.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace randlib_test
{

enum class Kind { Continuous, Discrete, Singular, BivariateContinuous, BivariateDiscrete };

/// Evaluates one recorded quantity, such as ("F", "0.5") or ("Mean", ""), the same way
/// tools/golden/dump_v1.cpp recorded it. Exceptions from the library propagate.
using Evaluator = std::function<double(const std::string &function, const std::string &arg)>;

struct Case
{
    std::string id;
    std::string className;
    Kind kind;
    std::function<Evaluator()> makeEvaluator;
    std::function<std::unique_ptr<ContinuousDistribution<double>>()> makeContinuous;
    std::function<std::unique_ptr<DiscreteDistribution<int>>()> makeDiscrete;
    std::function<std::unique_ptr<UnivariateDistribution<double>>()> makeSingular;
};

const std::vector<Case> &allCases();
const Case *findCase(const std::string &id);
std::vector<const Case *> casesOfKind(Kind kind);

}

#endif // RANDLIB_TESTS_REGISTRY_H
