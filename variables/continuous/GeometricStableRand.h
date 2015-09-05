#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "StableRand.h"
#include "ExponentialRand.h"
#include "LaplaceRand.h"

class RANDLIBSHARED_EXPORT GeometricStableRand : public StableRand
{
public:
    GeometricStableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual std::string name() override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

private:
    double variateForAlphaEqualOne() const;
    double variateForCommonAlpha() const;

public:
    void sample(QVector<double> &outputData);
};

#endif // GEOMETRICSTABLERAND_H
