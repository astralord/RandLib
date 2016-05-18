#ifndef NONCENTRALCHISQUARED_H
#define NONCENTRALCHISQUARED_H

#include "ContinuousDistribution.h"
#include "GammaRand.h"
#include "NormalRand.h"

/**
 * @brief The NoncentralChiSquared class
 */
class RANDLIBSHARED_EXPORT NoncentralChiSquared : public ContinuousDistribution
{
    double k, lambda;
    double sqrtLambda;

    ChiSquaredRand X;

public:
    explicit NoncentralChiSquared(int degree = 1, double noncentrality = 0);
    std::string name() override;

    void setParameters(int degree, double noncentrality);
    inline double getDegree() const { return k; }
    inline double getNoncentrality() const { return lambda; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // NONCENTRALCHISQUARED_H
