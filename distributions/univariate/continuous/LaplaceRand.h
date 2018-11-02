#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ExponentialRand.h"
#include "GeometricStableRand.h"

/**
 * @brief The AsymmetricLaplaceDistribution class <BR>
 * Abstract parent class for general Laplace and Asymmetric Laplace distributions
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT AsymmetricLaplaceDistribution : public GeneralGeometricStableDistribution<RealType>
{
public:
    AsymmetricLaplaceDistribution(double shift = 0, double scale = 1, double asymmetry = 1);

protected:
    void ChangeLocation();

public:
    void SetScale(double scale);

    inline double GetAsymmetry() const { return this->kappa; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    long double Entropy() const;

    void FitScale(const std::vector<RealType> &sample);
};


/**
 * @brief The CenteredAsymmetricLaplaceRand class <BR>
 * Centered Laplace distribution
 *
 * Belongs to exponential family
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT CenteredAsymmetricLaplaceRand : public AsymmetricLaplaceDistribution<RealType>,
                                                           public ExponentialFamily<RealType, DoublePair>
{
public:
    CenteredAsymmetricLaplaceRand(double scale = 1, double asymmetry = 1) :
        AsymmetricLaplaceDistribution<RealType>(0.0, scale, asymmetry) {}
};


/**
 * @brief The ShiftedAsymmetricLaplaceDistribution class <BR>
 * Abstract parent class for shifted Laplace and Asymmetric Laplace distributions
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT ShiftedAsymmetricLaplaceDistribution : public AsymmetricLaplaceDistribution<RealType>
{
public:
    ShiftedAsymmetricLaplaceDistribution(double shift = 0, double scale = 1, double asymmetry = 1) :
        AsymmetricLaplaceDistribution<RealType>(shift, scale, asymmetry) {}
    using GeneralGeometricStableDistribution<RealType>::SetShift;
    inline double GetShift() const { return this->m; }

    void FitShift(const std::vector<RealType> &sample);
protected:
    void FitShiftAndScale(const std::vector<RealType> &sample);
};


/**
 * @brief The AsymmetricLaplaceRand class <BR>
 * Asymmetric Laplace distribution
 *
 * Notation: X ~ Asymmetric-Laplace(m, γ, κ)
 *
 * Related distributions: <BR>
 * X = m + γ * (Y / κ - W * κ), where Y, W ~ Exp(1) <BR>
 * X - m ~ GS(2, β, γ, γ(1 - κ^2) / κ) with arbitrary β
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT AsymmetricLaplaceRand : public ShiftedAsymmetricLaplaceDistribution<RealType>
{
public:
    AsymmetricLaplaceRand(double shift = 0, double scale = 1, double asymmetry = 1) :
        ShiftedAsymmetricLaplaceDistribution<RealType>(shift, scale, asymmetry) {}
    String Name() const override;
    void SetAsymmetry(double asymmetry);
    static RealType StandardVariate(double asymmetry, RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

    using ShiftedAsymmetricLaplaceDistribution<RealType>::FitShiftAndScale;

private:
    DoublePair getOneSidedSums(const std::vector<RealType> &sample);
public:
    void FitAsymmetry(const std::vector<RealType> &sample);
};


/**
 * @brief The LaplaceRand class <BR>
 * Laplace distribution
 *
 * Notation: X ~ Laplace(m, γ)
 *
 * Related distributions: <BR>
 * X = m + γ * (Y - W), where Y, W ~ Exp(1) <BR>
 * X - m ~ GS(2, β, γ, 0) with arbitrary β
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT LaplaceRand : public ShiftedAsymmetricLaplaceDistribution<RealType>
{
public:
    LaplaceRand(double shift = 0, double scale = 1) :
        ShiftedAsymmetricLaplaceDistribution<RealType>(shift, scale, 1.0) {}
    String Name() const override;
    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);
    void Fit(const std::vector<RealType> &sample) { this->FitShiftAndScale(sample); }
};

#endif // LAPLACERAND_H
