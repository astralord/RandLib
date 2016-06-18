#ifndef CANTORRAND_H
#define CANTORRAND_H

#include "SingularDistribution.h"
#include "../discrete/BernoulliRand.h"

/**
 * @brief The CantorRand class
 */
class RANDLIBSHARED_EXPORT CantorRand : public SingularDistribution
{
    double generatorPrecision;
public:
    CantorRand();
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return 1; }

    void setGeneratorPrecision(double precision);
    inline double getGeneratorPrecision() const { return generatorPrecision; }

    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;
    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    double Likelihood(const std::vector<double> &sample) const override;
    double LogLikelihood(const std::vector<double> &sample) const override;
};

#endif // CANTORRAND_H
