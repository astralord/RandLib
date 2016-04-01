#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "LimitingDistribution.h"
#include "StableRand.h"

/**
 * @brief The GeometricStableRand class
 */
class RANDLIBSHARED_EXPORT GeometricStableRand : public LimitingDistribution
{
    StableRand Z;

public:
    GeometricStableRand(double exponent, double skewness, double scale = 1, double location = 0);
    std::string name() override;

    void setParameters(double exponent, double skewness);

    double f(double x) const override;
    double F(double x) const override;

private:
    double variateForAlphaEqualOne() const;
    double variateForAlphaEqualTwo() const;
    double variateForCommonAlpha() const;
public:
    double variate() const override;
    
    void sample(std::vector<double> &outputData) const override;
    
    std::complex<double> CF(double t) const override;

    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // GEOMETRICSTABLERAND_H
