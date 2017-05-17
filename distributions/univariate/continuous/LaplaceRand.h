#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ExponentialRand.h"
#include "GeometricStableRand.h"

/**
 * @brief The LaplaceRand class <BR>
 * Laplace distribution
 *
 * Notation: X ~ Laplace(m, γ, κ)
 *
 * Related distributions: <BR>
 * X = m + γ * (Y / κ - W * κ), where Y, W ~ Exp(1) <BR>
 * X - m ~ GeometricStable(2, β, γ, γ(1 - κ^2) / κ) with arbitrary β
 */
class RANDLIBSHARED_EXPORT LaplaceRand : public ShiftedGeometricStableDistribution
{
public:
    LaplaceRand(double shift = 0, double scale = 1, double asymmetry = 1);
    std::string Name() const override;

private:
    void ChangeLocation();
public:
    using ShiftedGeometricStableDistribution::SetShift;
    void SetScale(double scale);
    void SetAsymmetry(double asymmetry);

    inline double GetShift() const { return m; }
    inline double GetAsymmetry() const { return kappa; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;

    double Variate() const override;
    static double Variate(double location, double scale);
    static double Variate(double location, double scale, double asymmetry);
    void Sample(std::vector<double> &outputData) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;

    /// Maximum likelihood estimation
    /// One parameter
    void FitLocationMLE(const std::vector<double> &sample);
    void FitScaleMLE(const std::vector<double> &sample);
    void FitAsymmetryMLE(const std::vector<double> &sample);
    /// Two parameters
    void FitLocationAndScaleMLE(const std::vector<double> &sample);
    void FitLocationAndAsymmetryMLE(const std::vector<double> &sample);
    void FitScaleAndAsymmetryMLE(const std::vector<double> &sample);
    /// All parameters
    void FitLocationScaleAndAsymmetryMLE(const std::vector<double> &sample);

private:
    /**
     * @fn GetAsymmetryFromSkewness
     * @param sample
     * @return numeric solution of equation Skewness() = skewness
     */
    double GetAsymmetryFromSkewness(double skewness);
public:
    /// Method of moments
    /// One parameter
    void FitLocationMM(const std::vector<double> &sample);
    void FitScaleMM(const std::vector<double> &sample);
    void FitAsymmetryMM(const std::vector<double> &sample);
    /// Two parameters
    void FitLocationAndScaleMM(const std::vector<double> &sample);
    void FitLocationAndAsymmetryMM(const std::vector<double> &sample);
    void FitScaleAndAsymmetryMM(const std::vector<double> &sample);
    /// All parameters
    void FitLocationScaleAndAsymmetryMM(const std::vector<double> &sample);
};

#endif // LAPLACERAND_H
