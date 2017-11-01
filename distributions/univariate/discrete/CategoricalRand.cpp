#include "CategoricalRand.h"
#include "../continuous/UniformRand.h"

CategoricalRand::CategoricalRand(std::vector<double>&& probabilities)
{
    SetProbabilities(std::move(probabilities));
}

String CategoricalRand::Name() const
{
    String str = "Categorical(";
    for (int i = 0; i != K - 1; ++i)
        str += toStringWithPrecision(prob[i]) + ", ";
    return str + toStringWithPrecision(prob[K - 1]) + ")";
}

void CategoricalRand::SetProbabilities(std::vector<double> &&probabilities)
{
    if (probabilities.size() == 0 || std::accumulate(probabilities.begin(), probabilities.end(), 0.0) != 1.0)
        throw std::invalid_argument("Categorical distribution: probability parameters should sum to 1");
    else {
        prob = std::move(probabilities);
    }

    K = prob.size();
}

double CategoricalRand::P(const int & k) const
{
    return (k < 0 || k > K) ? 0.0 : prob[k];
}

double CategoricalRand::logP(const int & k) const
{
    return std::log(P(k));
}

double CategoricalRand::F(const int & k) const
{
    if (k < 0)
        return 0.0;
    if (k >= K)
        return 1.0;
    if (2 * k <= K) {
        double sum = 0.0;
        for (int i = 0; i <= k; ++i)
            sum += prob[i];
        return sum;
    }
    double sum = 1.0;
    for (int i = K - 1; i > k; --i)
        sum -= prob[i];
    return sum;
}

int CategoricalRand::Variate() const
{
    double U = UniformRand::StandardVariate();
    return quantileImpl(U);
}

double CategoricalRand::Mean() const
{
    double sum = 0.0;
    for (int i = 1; i != K; ++i)
        sum += i * prob[i];
    return sum;
}

double CategoricalRand::Variance() const
{
    double mean = 0.0, secMom = 0.0;
    for (int i = 1; i != K; ++i) {
        double aux = i * prob[i];
        mean += aux;
        secMom += i * aux;
    }
    return secMom - mean * mean;
}

int CategoricalRand::Mode() const
{
    auto maxProbIt = std::max_element(prob.begin(), prob.end());
    return std::distance(prob.begin(), maxProbIt);
}

int CategoricalRand::quantileImpl(double p) const
{
    double sum = 0.0;
    for (int i = 0; i != K; ++i) {
        sum += prob[i];
        if (sum >= p)
            return i;
    }
    return K;
}

int CategoricalRand::quantileImpl1m(double p) const
{
    double sum = prob[K - 1];
    for (int i = K - 2; i >= 0; --i) {
        sum += prob[i];
        if (sum >= p)
            return i;
    }
    return 0.0;
}

std::complex<double> CategoricalRand::CFImpl(double t) const
{
    double re = 0.0;
    double im = 0.0;
    for (int i = 0; i != K; ++i) {
        re += prob[i] * std::cos(t * i);
        im += prob[i] * std::sin(t * i);
    }
    return std::complex<double>(re, im);
}

