#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "StableRand.h"
#include "ExponentialRand.h"

class RANDLIBSHARED_EXPORT GeometricStableRand : public StableRand
{
public:
    GeometricStableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual std::string name() override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(QVector<double> &outputData);

    double E() const override { return (alpha > 1) ? mu : INFINITY; }
    double Var() const override { return (alpha == 2) ? 2 * sigma * sigma : INFINITY; }
};

#endif // GEOMETRICSTABLERAND_H
