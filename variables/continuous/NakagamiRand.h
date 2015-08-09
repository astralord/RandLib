#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

#include "ContinuousRand.h"
#include "GammaRand.h"

/**
 * @brief The NakagamiRand class
 */
class RANDLIBSHARED_EXPORT NakagamiRand : public ContinuousRand
{
    double m, w;
    double sigma; /// w / m
    double gammaMInv; /// 1.0 / gamma(m)
    double pdfCoef;

    GammaRand Y;

public:
    NakagamiRand(double shape = 0.5, double spread = 1);
    void setParameters(double shape, double spread);
    inline double getShape() const { return m; }
    inline double getSpread() const { return w; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() const override;

public:
    double E() const override { return std::tgamma(m + 0.5) * gammaMInv * std::sqrt(w / m); }
    double Var() const override {
        double res = std::tgamma(m + 0.5) * gammaMInv;
        res *= res;
        return w * (1 - res / m);
    }
};

#endif // NAKAGAMIRAND_H
