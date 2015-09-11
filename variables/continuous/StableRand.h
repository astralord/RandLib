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
protected:
    double alpha, beta, mu, sigma;

    // TODO: shouldn't storage them all
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
    virtual std::string name() override;

    void setParameters(double exponent, double skewness, double scale, double location);

    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }
    inline double getSigma() const { return sigma; }
    inline double getMu() const { return mu; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(QVector<double> &outputData);

    double E() const override { return (alpha > 1) ? mu : NAN; }
    double Var() const override { return (alpha == 2) ? 2 * sigma * sigma : INFINITY; }

private:
    /// variate
    double variateForCommonAlpha() const;
    double variateForAlphaEqualOne() const;

    /// pdf
    double pdfForCommonAlpha(double x) const;
    double integrandAuxForCommonAlpha(double theta, double xAdj, double xiAdj) const;
    double integrandForCommonAlpha(double theta, double xAdj, double xiAdj) const;

    double pdfForAlphaEqualOne(double x) const;
    double integrandAuxForAlphaEqualOne(double theta, double xAdj) const;
    double integrandForAlphaEqualOne(double theta, double xAdj) const;

    /// cdf
    // isn't written yet!
    double cdfForCommonAlpha(double x) const { return x; }
    double cdfForAlphaEqualOne(double x) const { return x; }
};


/**
 * @brief The HoltsmarkRand class
 */
class RANDLIBSHARED_EXPORT HoltsmarkRand : public StableRand
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableRand(1.5, 0.0, scale, location) {}
    virtual std::string name() override;
};

/**
 * @brief The LandauRand class
 */
class RANDLIBSHARED_EXPORT LandauRand : public StableRand
{
public:
    LandauRand(double scale = 1, double location = 0) : StableRand(1.0, 1.0, scale, location) {}
    virtual std::string name() override;
};

#endif // STABLERAND_H
