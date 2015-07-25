#ifndef STABLERAND_H
#define STABLERAND_H

#include <functional>
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"
#include "LevyRand.h"

/**
 * @brief The StableRand class
 */
class RANDLIBSHARED_EXPORT StableRand : public ContinuousRand
{
    double alpha, beta, mu, sigma;
    UniformRand U;
    ExponentialRand Exp;

    // TODO: shouldn't collect them all
    NormalRand N;
    CauchyRand C;
    LevyRand L;

    double B, S, alphaInv; /// coefficients for alpha != 1
    double logSigma; /// coefficients for alpha == 1
    double alpham1Inv, alpha_alpham1; /// 1 / (alpha - 1) and alpha / (alpha - 1)
    double pdfCoef;
    double zeta, xi, integrandCoef;

public:
    StableRand(double exponent, double skewness, double scale = 1, double location = 0);

    void setParameters(double exponent, double skewness, double scale, double location);

    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }
    inline double getSigma()  const{ return sigma; }
    inline double getMu() const { return mu; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() override;

    double E() const override { return (alpha > 1) ? mu : NAN; }
    double Var() const override { return (alpha == 2) ? 2 * sigma * sigma : INFINITY; }

private:
    /// value
    std::function<double ()> valuePtr;
    double valueForCommonAlpha();
    double valueForAlphaEqualOne();
    double valueNormal() { return N.variate(); }
    double valueCauchy() { return C.variate(); }
    double valueLevy() { return L.variate(); }
    double valueLevyNegative() { return -L.variate(); }

    /// pdf
    std::function<double (double)> pdfPtr;
    std::function<double (double)> integrandPtr; /// integrand for calculation of pdf

    double pdfForCommonAlpha(double x);
    double integrandAuxForCommonAlpha(double theta, double xAdj, double xiAdj) const;
    double rootFunctionForCommonAlpha(double theta, double xAdj, double xiAdj) const {
        return integrandAuxForCommonAlpha(theta, xAdj, xiAdj) - 1.0; }
    double integrandForCommonAlpha(double theta, double xAdj, double xiAdj) const;

    double pdfForAlphaEqualOne(double x);
    double integrandAuxForAlphaEqualOne(double theta, double xAdj) const;
    double rootFunctionForAlphaEqualOne(double theta, double xAdj) const {
        return integrandAuxForAlphaEqualOne(theta, xAdj) - 1.0; }
    double integrandForAlphaEqualOne(double theta, double xAdj) const;

    double pdfNormal(double x) { return N.f(x); }
    double pdfCauchy(double x) { return C.f(x); }
    double pdfLevy(double x) { return L.f(x); }
    double pdfLevyNegative(double x) { return L.f(-x); }

    /// cdf
    std::function<double (double)> cdfPtr;
    double cdfForCommonAlpha(double x) { return x; }
    double cdfForAlphaEqualOne(double x) { return x; }
    double cdfNormal(double x) { return N.F(x); }
    double cdfCauchy(double x) { return C.F(x); }
    double cdfLevy(double x) { return L.F(x); }
    double cdfLevyNegative(double x) { return 1.0 - L.F(-x); }
};

#endif // STABLERAND_H
