#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include <RandomVariable.h>

class RANDLIBSHARED_EXPORT UniformRand : public ContinuousRand
{
    double a, b;
    double c; /// c = 1 / (b - a)
    void swapBoundaries();

public:
    UniformRand(double minValue, double maxValue);
    void setBoundaries(double minValue, double maxValue);

    virtual double pdf (double x) const override;
    virtual double cdf(double x) const override;
    virtual double value() override;

    inline double M() const override { return .5 * (b + a); }
    inline double Var() const override { return (b - a) * (b - a) / 12; }

    inline double Median() const { return .5 * (b + a); }
    inline double Skewness() const { return 0; }
    inline double ExcessKurtosis() const { return -1.2; }

    inline double Entropy() const { return (b == a) ? -INFINITY : std::log(b - a); }

};

#endif // UNIFORMRAND_H
