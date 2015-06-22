#ifndef FISHERSNEDECORRAND_H
#define FISHERSNEDECORRAND_H

#include <RandomVariable.h>
#include "ChiSquaredRand.h"

class RANDLIBSHARED_EXPORT FisherSnedecorRand : public ContinuousRand
{
    int d1, d2;
    ChiSquaredRand U1, U2;

public:
    FisherSnedecorRand(int degree1, int degree2);

    void setDegrees(int degree1, int degree2);
    inline int getFirstDegree() const { return d1; }
    inline int getSecondDegree() const { return d2; }

    virtual double pdf (double x) const override;
    virtual double cdf(double x) const override;
    virtual double value() override;

    inline double M() const override { return (d2 > 2) ? d2 / (d2 - 2) : INFINITY /*or NAN*/; }
    inline double Var() const override {
        if (d2 <= 4)
            return INFINITY; /*or NAN*/
        double numen = 2 * d2 * d2 * (d1 + d2 - 2);
        double denom = d2 - 2;
        denom *= denom;
        denom *= d1 * (d2 - 4);
        return numen / denom;
    }
};

#endif // FISHERSNEDECORRAND_H
