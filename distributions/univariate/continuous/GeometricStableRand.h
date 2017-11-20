#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "StableRand.h"

/**
 * @brief The ShiftedGeometricStableDistribution class
 * Abstract class that unites in itself Asymmetric Laplace
 * and Geometric-Stable distributions
 */
template < typename RealType = long double >
class RANDLIBSHARED_EXPORT ShiftedGeometricStableDistribution : public ContinuousDistribution<RealType>
{
    StableRand<RealType> Z{};

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

    enum DISTRIBUTION_TYPE {
        LAPLACE, ///< α = 2, μ = 0
        ASYMMETRIC_LAPLACE, ///< α = 2, μ ≠ 0
        CAUCHY, ///< α = 1, β = 0
        LEVY, ///< α = 0.5, |β| = 1
        UNITY_EXPONENT, ///< α = 1, β ≠ 0
        ONEHALF_EXPONENT, ///< α = 0.5, |β| ≠ 1
        GENERAL ///< the rest
    };

    DISTRIBUTION_TYPE distributionType = LAPLACE; ///< type of distribution (Laplace by default)

public:
    inline double GetExponent() const { return alpha; }
    inline double GetSkewness() const { return beta; }
    inline double GetScale() const { return gamma; }
    inline double GetLogScale() const { return logGamma; }

    SUPPORT_TYPE SupportType() const override;
    RealType MinValue() const override;
    RealType MaxValue() const override;

protected:
    double pdfLaplace(double x) const;
    double logpdfLaplace(double x) const;
    double cdfLaplace(double x) const;
    double cdfLaplaceCompl(double x) const;

private:
    double pdfByLevy(double x) const;
    double pdfByCauchy(double x) const;

public:
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;

private:
    double variateForUnityExponent(double z) const;
    double variateForGeneralExponent(double z) const;
    double variateForOneHalfExponent(double z) const;
    double variateByCauchy(double z) const;
public:
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

protected:
    /**
     * @fn quantileLaplace
     * @param p input parameter in the interval (0, 1)
     * @return quantile for (Asymmetric) Laplace distribution
     */
    RealType quantileLaplace(double p) const;
    /**
     * @fn quantileLaplace1m
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for (Asymmetric) Laplace distribution
     */
    RealType quantileLaplace1m(double p) const;

    std::complex<double> CFImpl(double t) const override;
};


/**
 * @brief The GeometricStableRand class <BR>
 * Geometric-Stable distribution
 *
 * Notation: X ~ GS(α, β, γ, μ)
 *
 * If X ~ Laplace(m, γ, κ), then X - m ~ GS(2, β, γ, (1/κ - κ) * γ) with arbitrary β
 */
template < typename RealType = long double >
class RANDLIBSHARED_EXPORT GeometricStableRand : public ShiftedGeometricStableDistribution<RealType>
{
public:
    GeometricStableRand(double exponent, double skewness, double scale, double location);
    virtual ~GeometricStableRand() {}

    String Name() const override;
private:
    void ChangeAsymmetry();
public:
    void SetParameters(double exponent, double skewness);
    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return this->mu; }
};

#endif // GEOMETRICSTABLERAND_H
