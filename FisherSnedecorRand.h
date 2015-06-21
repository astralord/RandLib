#ifndef FISHERSNEDECORRAND_H
#define FISHERSNEDECORRAND_H

#include <RandomVariable.h>
#include <ChiSquaredRand.h>

class RANDLIBSHARED_EXPORT FisherSnedecorRand : public RandomVariable
{
    int d1, d2;
    ChiSquaredRand U1, U2;

public:
    FisherSnedecorRand(int degree1, int degree2);

    void setDegrees(int degree1, int degree2);
    int getFirstDegree() { return d1; }
    int getSecondDegree() { return d2; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    // TODO: do it
    inline double M() { return 0; }
    inline double Var() { return 0; }
};

#endif // FISHERSNEDECORRAND_H
