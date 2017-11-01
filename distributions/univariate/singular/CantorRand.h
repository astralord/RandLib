#ifndef CANTORRAND_H
#define CANTORRAND_H

#include "SingularDistribution.h"
#include "../discrete/BernoulliRand.h"

/**
 * @brief The CantorRand class <BR>
 * Cantor distribution
 *
 * Notation X ~ Cantor()
 */
class RANDLIBSHARED_EXPORT CantorRand : public SingularDistribution
{
    static constexpr int n = 30;
    static double table[n]; /// all powers of 1/3 from 1 to n
    static const bool dummy;
    static bool SetupTable();

public:
    CantorRand();
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return 1; }

    double F(const double & x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
};

#endif // CANTORRAND_H
