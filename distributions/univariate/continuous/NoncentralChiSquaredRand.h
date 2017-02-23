#ifndef NONCENTRALCHISQUAREDRAND_H
#define NONCENTRALCHISQUAREDRAND_H

#include "ContinuousDistribution.h"
#include "GammaRand.h"
#include "NormalRand.h"
#include "../discrete/PoissonRand.h"

/**
 * @brief The NoncentralChiSquaredRand class
 */
class RANDLIBSHARED_EXPORT NoncentralChiSquaredRand : public ContinuousDistribution
{
    double k, lambda;
    double halfK; /// 0.5 * k
    double sqrtLambda, logLambda;

    PoissonRand Y;

public:
    explicit NoncentralChiSquaredRand(double degree = 1, double noncentrality = 0);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double degree, double noncentrality);
    inline double GetDegree() const { return k; }
    inline double GetNoncentrality() const { return lambda; }

    double f(double x) const override;
    double logf(double x) const override;
    double F(double x) const override;
    double S(double x) const override;

private:
    double variateForDegreeEqualOne() const;

public:
    static double Variate(double degree, double noncentrality);
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // NONCENTRALCHISQUAREDRAND_H
