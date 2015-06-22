#ifndef PARETORAND_H
#define PARETORAND_H

#include <RandomVariable.h>
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousRand
{
    double xm, alpha;
    double alphaInv;
    double pdfCoef;
    UniformRand U;

public:
    ParetoRand(double shape, double scale);

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    double getShape() { return xm; }
    double getScale() { return alpha; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return (alpha > 1) ? alpha * xm / (alpha - 1) : INFINITY; }
    inline double Var() {
        if (alpha > 2)
        {
            double var = xm / (alpha - 1);
            var *= var;
            return alpha * var / (alpha - 2);
        }
        return (alpha > 1) ? INFINITY : NAN;
    }

    inline double Median() { return xm * std::pow(2.0, alphaInv); }
    inline double Mode() { return xm; }

    inline double Entropy() { return std::log(xm * alphaInv) + alphaInv + 1; }
};

#endif // PARETORAND_H
