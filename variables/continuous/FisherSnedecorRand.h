#ifndef FISHERSNEDECORRAND_H
#define FISHERSNEDECORRAND_H

#include "ContinuousRand.h"
#include "BetaPrimeRand.h"

/**
 * @brief The FisherSnedecorRand class
 */
class RANDLIBSHARED_EXPORT FisherSnedecorRand : public ContinuousRand
{
    int d1, d2;
    double pdfCoef;
    double a, d1_d2, c, d2_d1; /// constants for optimization

    BetaPrimeRand B;

public:
    FisherSnedecorRand(int degree1, int degree2);
    virtual std::string name() override;

    void setDegrees(int degree1, int degree2);
    void setFirstDegree(int degree1);
    void setSecondDegree(int degree2);
    inline int getFirstDegree() const { return d1; }
    inline int getSecondDegree() const { return d2; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(QVector<double> &outputData);

    double E() const override;
    double Var() const override;
};

#endif // FISHERSNEDECORRAND_H
