#ifndef CANTORRAND_H
#define CANTORRAND_H

#include "SingularRand.h"
#include "../discrete/BernoulliRand.h"

/**
 * @brief The CantorRand class
 */
class RANDLIBSHARED_EXPORT CantorRand : public SingularRand
{
    BernoulliRand B; /// for generator
    double generatorPrecision;
public:
    CantorRand();
    std::string name() override;

public:
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
};

#endif // CANTORRAND_H
