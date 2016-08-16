#ifndef CANTORRAND_H
#define CANTORRAND_H

#include "SingularDistribution.h"
#include "../discrete/BernoulliRand.h"

/**
 * @brief The CantorRand class
 * Cantor distribution
 *
 * Notation X ~ Cantor()
 */
class RANDLIBSHARED_EXPORT CantorRand : public SingularDistribution
{
    static constexpr int n = 30;
    static double table[n]; /// all powers of 1/3 from 1 to n
    static const bool dummy;
    static bool setupTable();

public:
    CantorRand();
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return 1; }

    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double QuantileImpl(double p) const override;
    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // CANTORRAND_H
