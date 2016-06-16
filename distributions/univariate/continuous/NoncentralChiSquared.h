#ifndef NONCENTRALCHISQUARED_H
#define NONCENTRALCHISQUARED_H

#include "ContinuousDistribution.h"
#include "GammaRand.h"
#include "NormalRand.h"
#include "../discrete/PoissonRand.h"

/**
 * @brief The NoncentralChiSquared class
 */
class RANDLIBSHARED_EXPORT NoncentralChiSquared : public ContinuousDistribution
{
    double k, lambda;
    double sqrtLambda;

    GammaRand X;
    PoissonRand Y;

public:
    explicit NoncentralChiSquared(double degree = 1, double noncentrality = 0);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void setParameters(double degree, double noncentrality);
    inline double getDegree() const { return k; }
    inline double getNoncentrality() const { return lambda; }

    double f(double x) const override;
    double F(double x) const override;

private:
    double variateForDegreeEqualOne() const;

public:
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // NONCENTRALCHISQUARED_H
