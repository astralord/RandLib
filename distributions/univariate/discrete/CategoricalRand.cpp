#include "CategoricalRand.h"
#include "../continuous/UniformRand.h"

CategoricalRand::CategoricalRand(std::vector<double>&& probabilities) :
    K(1), q(0.5)
{
    SetProbabilities(std::move(probabilities));
}

std::string CategoricalRand::Name() const
{
    std::string str = "Categorical(" + toStringWithPrecision(q) + ", ";
    for (double p : prob)
        str += toStringWithPrecision(p) + ", ";
    return str + toStringWithPrecision(K) + ")";
}

void CategoricalRand::SetProbabilities(std::vector<double> &&probabilities)
{
    long double sum = 0.0L;
    for (double var : probabilities) {
        sum += var;
    }
    if (probabilities.size() <= 0.0 || sum > 1.0) {
        /// we set default value,
        /// which is standard Bernoulli distribution
        prob = {0.5};
        sum = 0.5;
    }
    else {
        prob = std::move(probabilities);
    }

    K = prob.size();
    q = 1.0 - sum;
}

double CategoricalRand::P(int k) const
{
    if (k < 0 || k > K)
        return 0.0;
    return (k == 0) ? q : prob[k - 1];
}

double CategoricalRand::logP(int k) const
{
    return std::log(P(k));
}

double CategoricalRand::F(int k) const
{
    if (k < 0)
        return 0.0;
    if (k >= K)
        return 1.0;
    double sum = q;
    for (int i = 0; i != k; ++i) {
        sum += prob[i];
    }
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
    for (int i = 0; i != K; ++i)
        sum += (i + 1) * prob[i];
    return sum;
}

double CategoricalRand::Variance() const
{
    double mean = 0.0, secMom = 0.0;
    for (int i = 1; i <= K; ++i) {
        double aux = i * prob[i - 1];
        mean += aux;
        secMom += i * aux;
    }
    return secMom - mean * mean;
}

int CategoricalRand::Mode() const
{
    double maxProb = q;
    int mode = 0;
    for (int i = 0; i != K; ++i) {
        if (prob[i] > maxProb)
            mode = i + 1;
    }
    return mode;
}

double CategoricalRand::quantileImpl(double p) const
{
    double sum = q;
    for (int i = 0; i != K; ++i) {
        if (sum > p)
            return i;
        sum += prob[i];
    }
    return K;
}

double CategoricalRand::quantileImpl1m(double p) const
{
    double sum = prob[K - 1];
    for (int i = K - 2; i >= 0; --i) {
        if (sum > p)
            return i;
        sum += prob[i];
    }
    return 0.0;
}

std::complex<double> CategoricalRand::CFImpl(double t) const
{
    double re = q;
    double im = 0.0;
    for (int i = 0; i != K; ++i) {
        double ip1 = i + 1;
        re += prob[i] * std::cos(t * ip1);
        im += prob[i] * std::sin(t * ip1);
    }
    return std::complex<double>(re, im);
}

