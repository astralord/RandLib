#ifndef BETAPRIMERAND_H
#define BETAPRIMERAND_H

#include "BetaRand.h"

/**
 * @brief The BetaRand class
 */
class RANDLIBSHARED_EXPORT BetaPrimeRand : public BetaRand
{
public:
    BetaPrimeRand(double shape1, double shape2);
    virtual void setName() override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return (beta > 1) ? alpha / (beta - 1) : INFINITY; }
    double Var() const override {
        if (beta <= 2)
            return INFINITY;
        double betaAdj = beta - 1;
        double numerator = alpha * (alpha + betaAdj);
        double denominator = (betaAdj - 1) * betaAdj * betaAdj;
        return numerator / denominator;
    }
};

#endif // BETAPRIMERAND_H
