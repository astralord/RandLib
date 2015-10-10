#ifndef WEIBULLRAND_H
#define WEIBULLRAND_H

#include "ContinuousRand.h"

/**
 * @brief The WeibullRand class
 */
class RANDLIBSHARED_EXPORT WeibullRand : public ContinuousRand
{
    double lambda, k;
    double lambdaInv, kInv;

public:
    WeibullRand(double scale = 1, double shape = 1);
    std::string name() override;

    void setParameters(double scale, double shape);
    inline double getScale() const { return lambda; }
    inline double getShape() const { return k; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // WEIBULLRAND_H
