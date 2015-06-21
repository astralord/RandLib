#ifndef LEVYRAND_H
#define LEVYRAND_H

#include <RandomVariable.h>
#include <NormalRand.h>

/**
 * @brief The LevyRand class
 * X ~ Levy(mu, c)
 * X ~ Stable(0.5, 1, mu, c)
 */
class RANDLIBSHARED_EXPORT LevyRand : public RandomVariable
{
    double mu, c_2; /// c_2 = 0.5 * c
    double sqrtc_2pi; /// sqrtc_2pi = sqrt(0.5 * c / pi)
    NormalRand X;

public:
    LevyRand(double location, double scale);

    void setLocation(double location);
    void setScale(double scale);
    double getLocation() { return mu; }
    double getScale() { return c_2 + c_2; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return INFINITY; }
    inline double Var() { return INFINITY; }
};

/**
 * @brief The LevyNegativeRand class
 * -X ~ Levy(mu, c)
 *  X ~ Stable(0.5, -1, mu, c)
 */

class RANDLIBSHARED_EXPORT LevyNegativeRand : public LevyRand
{

public:
    LevyNegativeRand(double location, double scale) : LevyRand(-location, scale) {}

    void setLocation(double location) { LevyRand::setLocation(-location); }
    void setScale(double scale) { LevyRand::setScale(scale); }
    double getLocation() { return -LevyRand::getLocation(); }
    double getScale() { return LevyRand::getScale(); }

    virtual double pdf(double x) { return LevyRand::pdf(-x); }
    virtual double cdf(double x) { return 1.0 - LevyRand::pdf(-x); }
    virtual double value() { return -LevyRand::value(); }
};

#endif // LEVYRAND_H
