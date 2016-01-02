#ifndef ZETARAND_H
#define ZETARAND_H

#include "DiscreteRand.h"

/**
 * @brief The ZetaRand class
 *
 * P(X = k | n) = 1 / (k^s * zeta(s))
 */
class RANDLIBSHARED_EXPORT ZetaRand : public DiscreteRand
{
    double s, sm1;
    double zetaSInv; /// 1.0 / zeta(s)
    double b;
public:
    explicit ZetaRand(double exponent = 1.0);
    std::string name() override;

    void setExponent(double exponent);
    inline double getExponent() { return s; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;

    inline double getInverseZetaFunction() { return zetaSInv; }
};

#endif // ZETARAND_H
