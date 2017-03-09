#ifndef NONCENTRALTRAND_H
#define NONCENTRALTRAND_H

#include "ContinuousDistribution.h"
#include "StudentTRand.h"

/**
 * @brief The NoncentralTRand class
 * Notation: Noncentral-t(ν, μ)
 */
class RANDLIBSHARED_EXPORT NoncentralTRand : public ContinuousDistribution
{
    double nu, mu;
    double logNu;
    double PhiMu, PhimMu;
    double sqrt1p2oNu;

    struct nuCoefs {
        double qEpsCoef, q1mEpsCoef;
    } nuCoef, nup2Coef;

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
    DoublePair getIntegrationLimits(double x, double muAux, const nuCoefs &nuAuxCoef) const;
    double cdf(double x, double nuAux, double muAux, const nuCoefs &nuAuxCoef) const;
    double ccdf(double x, double nuAux, double muAux, const nuCoefs &nuAuxCoef) const;
    double g(double z, double x, double halfNuAux, double muAux, bool lower) const;
    double findMode(double x, double nuAux, double muAux, double A, double B) const;
    double lowerTail(const double & x, double nuAux, double muAux, const nuCoefs &nuAuxCoef, bool isCompl) const;
    double upperTail(const double & x, double nuAux, double muAux, const nuCoefs &nuAuxCoef, bool isCompl) const;

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // NONCENTRALTRAND_H
