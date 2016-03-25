#ifndef YULERAND_H
#define YULERAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/ParetoRand.h"

/**
 * @brief The YuleRand class
 * Yule distribution
 * X ~ Yule(\ro)
 */
class RANDLIBSHARED_EXPORT YuleRand : public DiscreteDistribution
{
    double ro;
    double gamma1pRo;
    
    ParetoRand X;
public:
    explicit YuleRand(double shape);
    std::string name() override;

    void setShape(double shape);
    inline double getShape() { return ro; }

    double P(int k) const override;
    double F(double x) const override;

    double variate() const override;
    static double variate(double shape);
    
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // YULERAND_H
