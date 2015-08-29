#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"
#include <functional>

/**
 * @brief The GeometricRand class
 */
class RANDLIBSHARED_EXPORT GeometricRand : public DiscreteRand<int>
{
    double p, q;

    static constexpr int tableSize = 16;
    // TODO: don't storage both variables (including tableSize)
    double table[tableSize];
    ExponentialRand W;

public:
    GeometricRand(double probability);
    virtual std::string name() override;

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData);

private:
    double variateForSmallP() const;
    double variateForLargeP() const;

public:
    double E() const override { return 1.0 / p; }
    double Var() const override { return (1 - p) / (p * p); }

    // TODO: add median and entropy
    static double constexpr Mode() { return 1; }
    inline double Skewness() { return (2 - p) / std::sqrt(1 - p); }
    inline double ExcessiveKurtosis() { return 1.0 / (p * (1 - p)) - 6; }
};

#endif // GEOMETRICRAND_H
