#ifndef FISHERSNEDECORRAND_H
#define FISHERSNEDECORRAND_H

#include "ContinuousRand.h"
#include "ChiSquaredRand.h"

/**
 * @brief The FisherSnedecorRand class
 */
class RANDLIBSHARED_EXPORT FisherSnedecorRand : public ContinuousRand
{
    int d1, d2;
    double gammaA, gammaB; /// gamma(.5 * d1) and gamma(.5 * d2)
    double pdfCoef;
    double a, b, c; /// constants for optimization
    ChiSquaredRand U1, U2;

public:
    FisherSnedecorRand(int degree1, int degree2);

    void setDegrees(int degree1, int degree2);
    void setFirstDegree(int degree1);
    void setSecondDegree(int degree2);
    inline int getFirstDegree() const { return d1; }
    inline int getSecondDegree() const { return d2; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return (d2 > 2) ? d2 / (d2 - 2) : INFINITY /*or NAN*/; }
    double Var() const override {
        if (d2 <= 4)
            return INFINITY; /*or NAN*/
        double numerator = 2 * d2 * d2 * (d1 + d2 - 2);
        double denominator = d2 - 2;
        denominator *= denominator;
        denominator *= d1 * (d2 - 4);
        return numerator / denominator;
    }
};

#endif // FISHERSNEDECORRAND_H
