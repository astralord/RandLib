#ifndef FISHERSNEDECORRAND_H
#define FISHERSNEDECORRAND_H

#include "BetaPrimeRand.h"

/**
 * @brief The FisherSnedecorRand class
 * F-distribution
 *
 * Notation: X ~ F(d1, d2)
 */
class RANDLIBSHARED_EXPORT FisherSnedecorRand : public ContinuousDistribution
{
    int d1, d2;
    double pdfCoef;
    double a, d1_d2, c, d2_d1; /// constants for optimization

    BetaPrimeRand B;

public:
    FisherSnedecorRand(int degree1, int degree2);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetDegrees(int degree1, int degree2);
    inline int GetFirstDegree() const { return d1; }
    inline int GetSecondDegree() const { return d2; }

    double f(double x) const override;
    double F(double x) const override;
    double S(double x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // FISHERSNEDECORRAND_H
