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
    explicit GeometricRand(double probability);
    virtual std::string name() override;

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double probability);

    void sample(QVector<double> &outputData);

private:
    double variateByExponential() const;
    double variateByTable() const;

public:
    double E() const override { return q / p; }
    double Var() const override { return q / (p * p); }

    double quantile(double F) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() const;
};

#endif // GEOMETRICRAND_H
