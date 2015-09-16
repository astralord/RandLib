#ifndef GUMBELRAND_H
#define GUMBELRAND_H

#include "ExponentialRand.h"

class RANDLIBSHARED_EXPORT GumbelRand : public ContinuousRand
{
    double mu, beta;
    double betaInv;
public:
    GumbelRand(double location, double scale);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double location, double scale);

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // GUMBELRAND_H
