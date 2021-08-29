#include "GeometricRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

template < typename IntType >
String GeometricRand<IntType>::Name() const
{
    return "Geometric(" + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template < typename IntType >
void GeometricRand<IntType>::SetProbability(double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("Geometric distribution: probability parameter should be in interval [0, 1]");
    this->SetParameters(1, probability);
}

template < typename IntType >
double GeometricRand<IntType>::P(const IntType &k) const
{
    return (k < 0) ? 0 : this->p * std::exp(k * this->log1mProb);
}

template < typename IntType >
double GeometricRand<IntType>::logP(const IntType &k) const
{
    return (k < 0) ? -INFINITY : this->logProb + k * this->log1mProb;
}

template < typename IntType >
double GeometricRand<IntType>::F(const IntType &k) const
{
    return (k < 0) ? 0 : -std::expm1l((k + 1) * this->log1mProb);
}

template < typename IntType >
double GeometricRand<IntType>::S(const IntType &k) const
{
    return (k < 0) ? 1 : std::exp((k + 1) * this->log1mProb);
}

template < typename IntType >
IntType GeometricRand<IntType>::Variate() const
{
    typename PascalRand<IntType>::GENERATOR_ID genId = this->GetIdOfUsedGenerator();
    if (genId == this->EXPONENTIAL)
        return this->variateGeometricThroughExponential();
    if (genId == this->TABLE)
        return this->variateGeometricByTable();
    /// unexpected return
    throw std::runtime_error("Geometric distribution: sampling failed");
}

template < typename IntType >
IntType GeometricRand<IntType>::Variate(double probability, RandGenerator &randGenerator)
{
    if (probability > 1.0 || probability < 0.0)
        throw std::invalid_argument("Geometric distribution: probability parameter should be in interval [0, 1], but it's equal to "
                                    + std::to_string(probability));

    /// here we use 0.05 instead of 0.08 because log(q) wasn't hashed
    if (probability < 0.05) {
        double rate = -std::log1pl(-probability);
        float X = ExponentialRand<float>::StandardVariate(randGenerator) / rate;
        return std::floor(X);
    }

    double U = UniformRand<double>::StandardVariate(randGenerator);
    int x = 0;
    double prod = probability, sum = prod, qprob = 1.0 - probability;
    while (U > sum) {
        prod *= qprob;
        sum += prod;
        ++x;
    }
    return x;
}

template < typename IntType >
void GeometricRand<IntType>::Sample(std::vector<IntType> &outputData) const
{
    typename PascalRand<IntType>::GENERATOR_ID genId = this->GetIdOfUsedGenerator();
    if (genId == this->EXPONENTIAL) {
        for (IntType &var : outputData)
            var = this->variateGeometricThroughExponential();
    }
    else if (genId == this->TABLE) {
        for (IntType &var : outputData)
            var = this->variateGeometricByTable();
    }
}

template < typename IntType >
IntType GeometricRand<IntType>::Median() const
{
    double median = -M_LN2 / this->log1mProb;
    double flooredMedian = std::floor(median);
    return (RandMath::areClose(median, flooredMedian, 1e-8)) ? flooredMedian - 1 : flooredMedian;
}

template < typename IntType >
long double GeometricRand<IntType>::Entropy() const
{
    double a = -this->q * this->log1mProb;
    double b = -this->p * this->logProb;
    return (a + b) / (M_LN2 * this->p);
}

template class GeometricRand<int>;
template class GeometricRand<long int>;
template class GeometricRand<long long int>;
