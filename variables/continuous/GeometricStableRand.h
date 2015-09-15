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

private:
    double variateForAlphaEqualOne() const;
    double variateForCommonAlpha() const;
public:
    double variate() const override;
    
    void sample(QVector<double> &outputData);
    
    std::complex<double> CF(double t) const override;
    
    double ExcessKurtosis() const override;
};

#endif // GEOMETRICSTABLERAND_H
