#include "registry.h"
#include "data.h"

#include <complex>
#include <stdexcept>

namespace randlib_test
{

namespace
{

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

std::string classOf(const std::string &id)
{
    return id.substr(0, id.find('('));
}

[[noreturn]] void unknown(const std::string &function)
{
    throw std::logic_error("test registry: unknown function " + function);
}

double complexPart(const std::string &function, std::complex<double> z)
{
    return function == "CF_re" ? z.real() : z.imag();
}

template <class D>
double univariateCommon(const D &d, const std::string &function, const std::string &arg, bool &handled)
{
    handled = true;
    if (function == "MinValue")
        return d.MinValue();
    if (function == "MaxValue")
        return d.MaxValue();
    if (function == "Mean")
        return d.Mean();
    if (function == "Variance")
        return d.Variance();
    if (function == "Skewness")
        return d.Skewness();
    if (function == "ExcessKurtosis")
        return d.ExcessKurtosis();
    if (function == "Median")
        return d.Median();
    if (function == "Quantile")
        return d.Quantile(parseNumber(arg));
    if (function == "Quantile1m")
        return d.Quantile1m(parseNumber(arg));
    if (function == "CF_re" || function == "CF_im")
        return complexPart(function, d.CF(parseNumber(arg)));
    handled = false;
    return 0.0;
}

Evaluator continuousEvaluator(std::shared_ptr<const ContinuousDistribution<double>> d)
{
    return [d](const std::string &function, const std::string &arg) {
        bool handled = false;
        double value = univariateCommon(*d, function, arg, handled);
        if (handled)
            return value;
        if (function == "Mode")
            return static_cast<double>(d->Mode());
        double x = arg.empty() ? 0.0 : parseNumber(arg);
        if (function == "f")
            return d->f(x);
        if (function == "logf")
            return d->logf(x);
        if (function == "F")
            return d->F(x);
        if (function == "S")
            return d->S(x);
        unknown(function);
    };
}

Evaluator discreteEvaluator(std::shared_ptr<const DiscreteDistribution<int>> d)
{
    return [d](const std::string &function, const std::string &arg) {
        bool handled = false;
        double value = univariateCommon(*d, function, arg, handled);
        if (handled)
            return value;
        if (function == "Mode")
            return static_cast<double>(d->Mode());
        int k = arg.empty() ? 0 : std::stoi(arg);
        if (function == "P")
            return d->P(k);
        if (function == "logP")
            return d->logP(k);
        if (function == "F")
            return d->F(k);
        if (function == "S")
            return d->S(k);
        unknown(function);
    };
}

Evaluator singularEvaluator(std::shared_ptr<const UnivariateDistribution<double>> d)
{
    return [d](const std::string &function, const std::string &arg) {
        bool handled = false;
        double value = univariateCommon(*d, function, arg, handled);
        if (handled)
            return value;
        double x = parseNumber(arg);
        if (function == "F")
            return d->F(x);
        if (function == "S")
            return d->S(x);
        unknown(function);
    };
}

template <class D>
double bivariateCommon(const D &d, const std::string &function, bool &handled)
{
    handled = true;
    if (function == "Mean1")
        return d.Mean().first;
    if (function == "Mean2")
        return d.Mean().second;
    if (function == "Cov11")
        return std::get<0>(d.Covariance());
    if (function == "Cov12")
        return std::get<1>(d.Covariance());
    if (function == "Cov22")
        return std::get<2>(d.Covariance());
    if (function == "Correlation")
        return d.Correlation();
    if (function == "Mode1")
        return d.Mode().first;
    if (function == "Mode2")
        return d.Mode().second;
    handled = false;
    return 0.0;
}

template <class D>
Evaluator bivariateContinuousEvaluator(std::shared_ptr<const D> d)
{
    return [d](const std::string &function, const std::string &arg) {
        bool handled = false;
        double value = bivariateCommon(*d, function, handled);
        if (handled)
            return value;
        auto semicolon = arg.find(';');
        DoublePair point(parseNumber(arg.substr(0, semicolon)), parseNumber(arg.substr(semicolon + 1)));
        if (function == "f")
            return d->f(point);
        if (function == "logf")
            return d->logf(point);
        unknown(function);
    };
}

template <class D>
Evaluator bivariateDiscreteEvaluator(std::shared_ptr<const D> d)
{
    return [d](const std::string &function, const std::string &arg) {
        bool handled = false;
        double value = bivariateCommon(*d, function, handled);
        if (handled)
            return value;
        auto semicolon = arg.find(';');
        IntPair point(std::stoi(arg.substr(0, semicolon)), std::stoi(arg.substr(semicolon + 1)));
        if (function == "P")
            return d->P(point);
        if (function == "logP")
            return d->logP(point);
        unknown(function);
    };
}

std::vector<Case> buildCases()
{
    std::vector<Case> cases;

#define RL_PLAIN(X) X
#define RL_CONT(T, ...)                                                                                                \
    cases.push_back({caseId(#T, #__VA_ARGS__), classOf(caseId(#T, #__VA_ARGS__)), Kind::Continuous,                    \
                     [] { return continuousEvaluator(std::make_shared<const T>(__VA_ARGS__)); },                       \
                     [] { return std::unique_ptr<ContinuousDistribution<double>>(new T(__VA_ARGS__)); }, nullptr,      \
                     nullptr});
#define RL_DISC(T, ...)                                                                                                \
    cases.push_back({caseId(#T, #__VA_ARGS__), classOf(caseId(#T, #__VA_ARGS__)), Kind::Discrete,                      \
                     [] { return discreteEvaluator(std::make_shared<const T>(__VA_ARGS__)); }, nullptr,                \
                     [] { return std::unique_ptr<DiscreteDistribution<int>>(new T(__VA_ARGS__)); }, nullptr});
#define RL_SING(T, ...)                                                                                                \
    cases.push_back({caseId(#T, #__VA_ARGS__), classOf(caseId(#T, #__VA_ARGS__)), Kind::Singular,                      \
                     [] { return singularEvaluator(std::make_shared<const T>(__VA_ARGS__)); }, nullptr, nullptr,       \
                     [] { return std::unique_ptr<UnivariateDistribution<double>>(new T(__VA_ARGS__)); }});
#define RL_BIVC(T, ...)                                                                                                \
    cases.push_back({caseId(#T, #__VA_ARGS__), classOf(caseId(#T, #__VA_ARGS__)), Kind::BivariateContinuous,           \
                     [] { return bivariateContinuousEvaluator(std::make_shared<const T>(__VA_ARGS__)); }, nullptr,     \
                     nullptr, nullptr});
#define RL_BIVD(T, ...)                                                                                                \
    cases.push_back({caseId(#T, #__VA_ARGS__), classOf(caseId(#T, #__VA_ARGS__)), Kind::BivariateDiscrete,             \
                     [] { return bivariateDiscreteEvaluator(std::make_shared<const T>(__VA_ARGS__)); }, nullptr,       \
                     nullptr, nullptr});

#include "cases.inc"

#undef RL_PLAIN
#undef RL_CONT
#undef RL_DISC
#undef RL_SING
#undef RL_BIVC
#undef RL_BIVD
    return cases;
}

} // namespace

const std::vector<Case> &allCases()
{
    static const std::vector<Case> cases = buildCases();
    return cases;
}

const Case *findCase(const std::string &id)
{
    for (const Case &c : allCases()) {
        if (c.id == id)
            return &c;
    }
    return nullptr;
}

std::vector<const Case *> casesOfKind(Kind kind)
{
    std::vector<const Case *> result;
    for (const Case &c : allCases()) {
        if (c.kind == kind)
            result.push_back(&c);
    }
    return result;
}

} // namespace randlib_test
