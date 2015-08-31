#ifndef PARETORAND_H
#define PARETORAND_H

#include "ContinuousRand.h"
#include "ExponentialRand.h"


/**
 * @brief The ParetoRand class
 */
class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousRand
{
    double xm, alpha;
    double alphaInv;
    double pdfCoef;

public:
    ParetoRand(double shape, double scale);
    virtual std::string name() override;

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    inline double getShape() const { return xm; }
    inline double getScale() const { return alpha; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double shape, double scale);

private:
    static double variateForAlphaOne();
    static double variateForAlphaTwo();
    static double variateForCommonAlpha(double scale);

public:
    void sample(QVector<double> &outputData);

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
