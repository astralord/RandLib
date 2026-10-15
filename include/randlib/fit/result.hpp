#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace randlib
{

enum class Method
{
    MLE,
    MethodOfMoments,
    UMVU,
    Bayes,
    Quantile
};

enum class Criterion
{
    AIC,
    BIC
};

struct FitOptions
{
    Method method = Method::MLE;
    std::map<std::string, double, std::less<>> fixed;
    double tolerance = 1e-8;
    int max_iterations = 20000;
};

template <class Distribution>
struct FitResult
{
    explicit FitResult(Distribution fittedDistribution) : distribution(std::move(fittedDistribution)) {}

    Distribution distribution;
    double log_likelihood = 0;
    double aic = 0;
    double bic = 0;
    std::vector<std::string> parameter_names;
    std::vector<double> parameters;
    std::vector<double> std_errors;
    std::vector<std::vector<double>> covariance;
    std::vector<bool> fixed;
    bool converged = false;
    int iterations = 0;
    Method method = Method::MLE;
};

} // namespace randlib
