#ifndef YULERAND_H
#define YULERAND_H

#include "DiscreteRand.h"
#include "GeometricRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/ParetoRand.h"

/**
 * @brief The YuleRand class
 */
class RANDLIBSHARED_EXPORT YuleRand : public DiscreteRand
{
    double ro;
    double gamma1pRo;
public:
    explicit YuleRand(double shape);
    virtual std::string name() override;

    void setShape(double shape);
    inline double getShape() { return ro; }

    double P(int k) const override;
    double F(double x) const override;

    double variate() const override;
    static double variate(double shape);

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // YULERAND_H
