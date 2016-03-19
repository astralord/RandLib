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
class RANDLIBSHARED_EXPORT StableRand : public ContinuousDistribution
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

    // TODO: find an appropriate value
    static constexpr double almostPI_2 = 1.555; /// ~ 0.99 * (pi / 2)

public:
    StableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~StableRand() {}
    std::string name() override;

    void setParameters(double exponent, double skewness, double scale, double location);

    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }
    inline double getSigma() const { return sigma; }
    inline double getMu() const { return mu; }
    
private:
    /// pdf
    double pdfForCommonAlpha(double x) const;
    double integrandAuxForCommonAlpha(double theta, double xAdj, double xiAdj) const;
    double integrandForCommonAlpha(double theta, double xAdj, double xiAdj) const;

    double pdfForAlphaEqualOne(double x) const;
    double integrandAuxForAlphaEqualOne(double theta, double xAdj) const;
    double integrandForAlphaEqualOne(double theta, double xAdj) const;
public:    
    double f(double x) const override;
    
private:
    /// cdf
    double cdfForCommonAlpha(double x) const;
    double cdfForAlphaEqualOne(double x) const;
public:
    double F(double x) const override;
    
private:
    /// variate
    double variateForCommonAlpha() const;
    double variateForAlphaEqualOne() const;
public:
    double variate() const override;
    void sample(QVector<double> &outputData) const override;

protected:
    std::complex<double> psi(double t) const;
public:
    std::complex<double> CF(double t) const override;

    double Mean() const override;
    double Variance() const override;

    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The HoltsmarkRand class
 * X ~ Stable(1.5, 0.0, scale, location)
 */
class RANDLIBSHARED_EXPORT HoltsmarkRand : public StableRand
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableRand(1.5, 0.0, scale, location) {}
    std::string name() override;
};

/**
 * @brief The LandauRand class
 * X ~ Stable(1.0, 1.0, scale, location)
 */
class RANDLIBSHARED_EXPORT LandauRand : public StableRand
{
public:
    LandauRand(double scale = 1, double location = 0) : StableRand(1.0, 1.0, scale, location) {}
    std::string name() override;
};

#endif // STABLERAND_H
