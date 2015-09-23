#ifndef ZETARAND_H
#define ZETARAND_H

#include "DiscreteRand.h"

class RANDLIBSHARED_EXPORT ZetaRand : public DiscreteRand
{
    double s;
    double zetaSInv; /// 1.0 / zeta(s)
    double b;
public:
    explicit ZetaRand(double exponent);
    std::string name() override;

    void setExponent(double exponent);
    inline double getExponent() { return s; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
};

#endif // ZETARAND_H
