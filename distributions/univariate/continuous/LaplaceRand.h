#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ExponentialRand.h"
#include "GeometricStableRand.h"

/**
 * @brief The AsymmetricLaplaceDistribution class <BR>
 * Abstract parent class for Laplace and Asymmetric Laplace distributions
 */
class RANDLIBSHARED_EXPORT AsymmetricLaplaceDistribution : public ShiftedGeometricStableDistribution
{
public:
    AsymmetricLaplaceDistribution(double shift = 0, double scale = 1, double asymmetry = 1);

protected:
    void ChangeLocation();

public:
    using ShiftedGeometricStableDistribution::SetShift;
    void SetScale(double scale);

    inline double GetShift() const { return m; }
    inline double GetAsymmetry() const { return kappa; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;

    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;

    void FitShift(const std::vector<double> &sample);
    void FitScale(const std::vector<double> &sample);

protected:
    void FitShiftAndScale(const std::vector<double> &sample);
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
class RANDLIBSHARED_EXPORT AsymmetricLaplaceRand : public AsymmetricLaplaceDistribution
{
public:
    AsymmetricLaplaceRand(double shift = 0, double scale = 1, double asymmetry = 1) : AsymmetricLaplaceDistribution(shift, scale, asymmetry) {}
    String Name() const override;
    void SetAsymmetry(double asymmetry);
    static double StandardVariate(double asymmetry);

    using AsymmetricLaplaceDistribution::FitShiftAndScale;

    void FitAsymmetry(const std::vector<double> &sample);
    void FitShiftAndAsymmetry(const std::vector<double> &sample);
    void FitScaleAndAsymmetry(const std::vector<double> &sample);
    void Fit(const std::vector<double> &sample);
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
class RANDLIBSHARED_EXPORT LaplaceRand : public AsymmetricLaplaceDistribution
{
public:
    LaplaceRand(double shift = 0, double scale = 1) : AsymmetricLaplaceDistribution(shift, scale, 1.0) {}
    String Name() const override;
    static double StandardVariate();
    void Fit(const std::vector<double> &sample) { FitShiftAndScale(sample); }
};

#endif // LAPLACERAND_H
