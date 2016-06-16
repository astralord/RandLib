#ifndef FISHERSNEDECORRAND_H
#define FISHERSNEDECORRAND_H

#include "BetaPrimeRand.h"

/**
 * @brief The FisherSnedecorRand class
 * F-distribution
 * X ~ F(d_1, d_2)
 */
class RANDLIBSHARED_EXPORT FisherSnedecorRand : public ContinuousDistribution
{
    int d1, d2;
    double pdfCoef;
    double a, d1_d2, c, d2_d1; /// constants for optimization

    BetaPrimeRand B;

public:
    FisherSnedecorRand(int degree1, int degree2);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return SEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void setDegrees(int degree1, int degree2);
    void setFirstDegree(int degree1);
    void setSecondDegree(int degree2);
    inline int getFirstDegree() const { return d1; }
    inline int getSecondDegree() const { return d2; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;
    
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // FISHERSNEDECORRAND_H
