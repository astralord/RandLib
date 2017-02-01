#ifndef NONCENTRALTRAND_H
#define NONCENTRALTRAND_H

#include "ContinuousDistribution.h"
#include "StudentTRand.h"

/**
 * @brief The NoncentralTRand class
 */
class RANDLIBSHARED_EXPORT NoncentralTRand : public ContinuousDistribution
{
    double nu, mu;
    double logNu;
    double cdfCoef;

    StudentTRand T;

public:
    explicit NoncentralTRand(double degree = 1, double noncentrality = 0);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double degree, double noncentrality);
    inline double GetDegree() const { return nu; }
    inline double GetNoncentrality() const { return mu; }

private:
    double Faux(double x, double nuAux, double muAux) const;

public:
    double f(double x) const override;
    double F(double x) const override;

    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // NONCENTRALTRAND_H
