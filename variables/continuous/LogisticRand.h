#ifndef LOGISTICRAND_H
#define LOGISTICRAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT LogisticRand : public ContinuousRand
{
    double mu, s;
    UniformRand U;

public:
    LogisticRand(double location = 0, double scale = 1);

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() { return mu; }
    inline double getScale() { return s; }

    virtual double f(double x) const;
    virtual double F(double x) const;
    virtual double value();

    double E() const { return mu; }
    double Var() const {
        double sPi = s * M_PI;
        return sPi * sPi / 3;
    }

    inline double Median() const { return mu; }
    inline double Mode() const { return mu; }
    static constexpr double Skewness() { return 0; }
    static constexpr double ExcessKurtosis() { return 1.2; }
};

#endif // LOGISTICRAND_H
