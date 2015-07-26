#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include "ContinuousRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The UniformRand class
 * Continuous uniform distribution: X ~ U(a, b)
 */
class RANDLIBSHARED_EXPORT UniformRand : public ContinuousRand
{
    double a, b;
    double c; /// 1 / (b - a)
    double delta; /// (b - a) / RAND_MAX
    void swapBoundaries();

public:
    UniformRand(double minValue = 0, double maxValue = 1);

    void setBoundaries(double minValue, double maxValue);    
    inline double getA() { return a; }
    inline double getB() { return b; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() override;

    double E() const override { return .5 * (b + a); }
    double Var() const override { return (b - a) * (b - a) / 12; }

    inline double Median() const { return .5 * (b + a); }
    static constexpr double Skewness() { return 0; }
    static constexpr double ExcessKurtosis() { return -1.2; }

    inline double Entropy() const { return (b == a) ? -INFINITY : std::log(b - a); }

    bool fitToData(const QVector<double> &sample);
};

#endif // UNIFORMRAND_H
