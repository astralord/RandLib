#ifndef PARETORAND_H
#define PARETORAND_H

#include "ContinuousRand.h"
#include "ExponentialRand.h"

// TODO: maybe E/Erlang is faster for generator
class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousRand
{
    double xm, alpha;
    double alphaInv;
    double pdfCoef;

public:
    ParetoRand(double shape, double scale);
    virtual void setName() override;

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    inline double getShape() const { return xm; }
    inline double getScale() const { return alpha; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() const override;

    double E() const override { return (alpha > 1) ? alpha * xm / (alpha - 1) : INFINITY; }
    double Var() const override {
        if (alpha > 2)
        {
            double var = xm / (alpha - 1);
            var *= var;
            return alpha * var / (alpha - 2);
        }
        return (alpha > 1) ? INFINITY : NAN;
    }

    inline double Median() const { return xm * std::pow(2.0, alphaInv); }
    inline double Mode() const { return xm; }

    inline double Entropy() const { return std::log(xm * alphaInv) + alphaInv + 1; }

    bool fitToData(const QVector<double> &sample);
};

#endif // PARETORAND_H
