#ifndef YULERAND_H
#define YULERAND_H

#include "DiscreteRand.h"
#include "GeometricRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/ParetoRand.h"

/**
 * @brief The YuleRand class
 */
class RANDLIBSHARED_EXPORT YuleRand : public DiscreteRand<int>
{
    double ro;
    double gamma1pRo;
public:
    YuleRand(double shape);
    virtual std::string name() override;

    void setShape(double shape);
    inline double getShape() { return ro; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double shape);

public:
    double E() const override { return (ro <= 1) ? INFINITY : ro / (ro - 1); }
    double Var() const override {
        if (ro <= 2)
            return INFINITY;
        double aux = ro / (ro - 1);
        return aux * aux / (ro - 2);
    }

    static constexpr double Mode() { return 1.0; }
};

#endif // YULERAND_H
