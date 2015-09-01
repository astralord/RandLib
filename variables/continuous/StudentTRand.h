#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "ContinuousRand.h"
#include "ChiSquaredRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"

/**
 * @brief The StudentTRand class
 */
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousRand
{
    int v;
    ChiSquaredRand Y;
    double pdfCoef;
public:
    StudentTRand(int degree);
    virtual std::string name() override;

    void setDegree(int degree);
    inline int getDegree() const { return v; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData);

    double E() const override { return (v > 1) ? 0 : NAN; }
    double Var() const override {
        if (v > 2)
            return v / (v - 2);
        return (v > 1) ? INFINITY : NAN;
    }
};

#endif // STUDENTTRAND_H
