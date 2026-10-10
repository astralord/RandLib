#include "CategoricalRand.h"
#include "../continuous/UniformRand.h"

template < typename IntType >
CategoricalRand<IntType>::CategoricalRand(std::vector<double>&& probabilities)
{
    SetProbabilities(std::move(probabilities));
}

template < typename IntType >
String CategoricalRand<IntType>::Name() const
{
    String str = "Categorical(";
    for (int i = 0; i != K - 1; ++i)
        str += this->toStringWithPrecision(prob[i]) + ", ";
    return str + this->toStringWithPrecision(prob[K - 1]) + ")";
}

template < typename IntType >
void CategoricalRand<IntType>::SetProbabilities(std::vector<double> &&probabilities)
{
    if (probabilities.empty())
        throw std::invalid_argument("Categorical distribution: there should be at least one probability");
    for (double probability : probabilities) {
        if (!(probability >= 0.0))
            throw std::invalid_argument("Categorical distribution: probabilities should be non-negative");
    }
    double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    if (!(std::fabs(sum - 1.0) <= 1e-12 * probabilities.size()))
        throw std::invalid_argument("Categorical distribution: probability parameters should sum to 1");
    prob = std::move(probabilities);
    K = prob.size();
}

template < typename IntType >
double CategoricalRand<IntType>::P(const IntType & k) const
{
    return (k < 0 || k >= K) ? 0.0 : prob[k];
}

template < typename IntType >
double CategoricalRand<IntType>::logP(const IntType & k) const
{
    return std::log(P(k));
}

template < typename IntType >
double CategoricalRand<IntType>::F(const IntType & k) const
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

template < typename IntType >
IntType CategoricalRand<IntType>::Variate() const
{
    double U = UniformRand<double>::StandardVariate(this->localRandGenerator);
    return quantileImpl(U);
}

template < typename IntType >
long double CategoricalRand<IntType>::Mean() const
{
    long double sum = 0.0;
    for (int i = 1; i != K; ++i)
        sum += i * prob[i];
    return sum;
}

template < typename IntType >
long double CategoricalRand<IntType>::Variance() const
{
    long double mean = 0.0, secMom = 0.0;
    for (int i = 1; i != K; ++i) {
        long double aux = i * prob[i];
        mean += aux;
        secMom += i * aux;
    }
    return secMom - mean * mean;
}

template < typename IntType >
IntType CategoricalRand<IntType>::Mode() const
{
    auto maxProbIt = std::max_element(prob.begin(), prob.end());
    return std::distance(prob.begin(), maxProbIt);
}

template < typename IntType >
IntType CategoricalRand<IntType>::quantileImpl(double p) const
{
    /// smallest k with F(k) ≥ p
    double sum = 0.0;
    for (IntType i = 0; i != K; ++i) {
        sum += prob[i];
        if (sum >= p)
            return i;
    }
    return K - 1;
}

template < typename IntType >
IntType CategoricalRand<IntType>::quantileImpl1m(double p) const
{
    /// smallest k with S(k) ≤ p, where S(k - 1) = S(k) + P(k)
    double tail = 0.0;
    for (IntType k = K - 1; k > 0; --k) {
        tail += prob[k];
        if (tail > p)
            return k;
    }
    return 0;
}

template < typename IntType >
std::complex<double> CategoricalRand<IntType>::CFImpl(double t) const
{
    double re = 0.0;
    double im = 0.0;
    for (int i = 0; i != K; ++i) {
        re += prob[i] * std::cos(t * i);
        im += prob[i] * std::sin(t * i);
    }
    return std::complex<double>(re, im);
}


template class CategoricalRand<int>;
template class CategoricalRand<long int>;
template class CategoricalRand<long long int>;
