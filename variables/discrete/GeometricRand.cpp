#include "GeometricRand.h"

GeometricRand::GeometricRand(double probability)
{
    setProbability(probability);
}

std::string GeometricRand::name()
{
    return "Geometric(" + toStringWithPrecision(getProbability()) + ")";
}

void GeometricRand::setProbability(double probability)
{
    p = std::min(probability, 1.0);
    if (p <= 0)
        p = MIN_POSITIVE;
    q = 1.0 - p;

    /// we use two different generators for two different cases
    /// if p < 0.2 then the tail is too heavy
    /// (probability to be in main body is less then 0.977482)
    /// thus we always use low bound of variate exponential distribution
    /// otherwise we choose table method
    if (p < 0.2)
    {
        W.setRate(-std::log(q));
    }
    else
    {
        table[0] = p;
        double prod = p;
        for (int i = 1; i < tableSize; ++i)
        {
            prod *= q;
            table[i] = table[i - 1] + prod;
        }
    }
}

double GeometricRand::P(int k) const
{
    return (k < 0) ? 0 : p * std::pow(q, k);
}

double GeometricRand::F(double x) const
{
    return (x < 0) ? 0 : 1 - std::pow(q, std::floor(x) + 1);
}

double GeometricRand::variate() const
{
    if (p < 0.2)
        return variateByExponential();
    return variateByTable();
}

double GeometricRand::variate(double probability)
{
    if (probability < 0.2)
        return std::floor(ExponentialRand::variate(-std::log(1 - probability)));

    double U = UniformRand::standardVariate();
    int x = 0;
    double prod = probability, sum = prod, q = 1 - probability;
    while (U > sum) {
        prod *= q;
        sum += prod;
        ++x;
    }
    return x;
}

void GeometricRand::sample(QVector<double> &outputData)
{
    if (p < 0.2) {
        for (double &var : outputData)
            var = variateByExponential();
    }
    else {
        for (double &var : outputData)
            var = variateByTable();
    }
}

double GeometricRand::variateByExponential() const
{
    return std::floor(W.variate());
}

double GeometricRand::variateByTable() const
{
    double U = UniformRand::standardVariate();
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateByTable();
    /// handle the main body
    int x = 0;
    while (U > table[x])
        ++x;
    return x;
}

double GeometricRand::Quantile(double F) const
{
    if (F < 0 || F > 1)
        return NAN;
    if (F == 0)
        return -INFINITY;
    if (F == 1)
        return INFINITY;
    return std::log1p(F) / std::log(1 - p);
}

double GeometricRand::Median() const
{
    return std::floor(-M_LN2 / std::log(q));
}

double GeometricRand::Mode() const
{
    return 0.0;
}

double GeometricRand::Skewness() const
{
    return (2 - p) / std::sqrt(q);
}

double GeometricRand::ExcessKurtosis() const
{
    return p * p / q + 6;
}

double GeometricRand::Entropy() const
{
    double a = -q * std::log(q);
    double b = -p * std::log(p);
    return (a + b) / (M_LN2 * p);
}

