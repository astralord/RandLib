#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

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
    std::string name() override;

    void setParameters(double shape, double spread);
    inline double getShape() const { return m; }
    inline double getSpread() const { return w; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

public:
    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
};

#endif // NAKAGAMIRAND_H
