#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "StableRand.h"

/**
 * @brief The ShiftedGeometricStableDistribution class
 * Abstract class that unites in itself Asymmetric Laplace
 * and Geometric-Stable distributions
 */
class RANDLIBSHARED_EXPORT ShiftedGeometricStableDistribution : public ContinuousDistribution
{
    StableRand Z{};

protected:
    double alpha = 2; ///< characteristic exponent α
    double alphaInv = 0.5; /// 1/α
    double beta = 0; ///< skewness β
    double mu = 0; ///< location μ
    double m = 0; ///< shift m
    double gamma = M_SQRT2; ///< scale γ
    double logGamma = 0.5 * M_LN2; ///< log(γ)

    /// parameters for α = 2
    double kappa = 1; ///< asymmetry coefficient κ
    double kappaInv = 1; ///< 1 / κ
    double kappaSq = 1; ///< κ^2
    double log1pKappaSq = M_LN2; ///< log(1 + κ^2)
    double pdfCoef = M_LN2; ///< log(γ * (κ + 1 / κ))
    double cdfCoef = -M_LN2; ///< 2 * log(κ) - log(1 + κ^2)

    ShiftedGeometricStableDistribution(double exponent, double skewness, double scale = 1.0, double location = 0.0, double shift = 0.0);
    virtual ~ShiftedGeometricStableDistribution() {}

    void SetParameters(double exponent, double skewness, double scale = 1.0, double location = 0.0, double shift = 0.);
    void SetLocation(double location);
    void SetShift(double shift);
    void SetScale(double scale);
    void SetAsymmetry(double asymmetry);

public:
    inline double GetExponent() const { return alpha; }
    inline double GetSkewness() const { return beta; }
    inline double GetScale() const { return gamma; }
    inline double GetLogScale() const { return logGamma; }

    SUPPORT_TYPE SupportType() const override;
    double MinValue() const override;
    double MaxValue() const override;

protected:
    double pdfLaplace(double x) const;
    double logpdfLaplace(double x) const;
    double cdfLaplace(double x) const;
    double cdfLaplaceCompl(double x) const;

private:
    double pdfByLevy(double x) const;
    double pdfByCauchy(double x) const;

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;

private:
    double variateForUnityExponent() const;
    double variateForCommonExponent() const;
    double variateByLevy(bool positive) const;
    double variateByCauchy() const;
public:
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};


/**
 * @brief The GeometricStableRand class <BR>
 * Geometric-Stable distribution
 *
 * Notation: X ~ Geometric-Stable(α, β, γ, μ)
 *
 * If X ~ Laplace(m, γ, κ), then X - m ~ Geometric-Stable(2, β, γ, (1/κ - κ) * γ) with arbitrary β
 */
class RANDLIBSHARED_EXPORT GeometricStableRand : public ShiftedGeometricStableDistribution
{
public:
    GeometricStableRand(double exponent, double skewness, double scale, double location) : ShiftedGeometricStableDistribution(exponent, skewness, scale, location) {}
    virtual ~GeometricStableRand() {}

    std::string Name() const override;
private:
    void ChangeAsymmetry();
public:
    void SetParameters(double exponent, double skewness);
    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return mu; }
};

#endif // GEOMETRICSTABLERAND_H
