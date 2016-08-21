#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ExponentialRand.h"
#include "GeometricStableRand.h"

/**
 * @brief The LaplaceRand class
 * Laplace distribution
 *
 * Notation: X ~ Laplace(m, σ, k)
 *
 * Related distributions:
 * X = m + σ (Y / k - W * k), where Y, W ~ Exp(1)
 * If X ~ Laplace(m, σ, k), then X - m ~ GeometricStable(2, 0, σ, σ(1 - k^2) / k)
 */
class RANDLIBSHARED_EXPORT LaplaceRand : public GeometricStableRand
{
    double m;

public:
    LaplaceRand(double shift = 0, double scale = 1, double asymmetry = 1);
    std::string Name() const override;

    void SetShift(double shift);
    void SetAsymmetry(double asymmetry);
    inline double GetShift() const { return m; }
    inline double GetAsymmetry() const { return k; }

    double f(double x) const override;
    double F(double x) const override;

    double Variate() const override;
    static double Variate(double location, double scale);
    static double Variate(double location, double scale, double asymmetry);
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Median() const override;
    double Mode() const override;

    std::complex<double> CF(double t) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double Entropy() const;

    /// Maximum likelihood estimation
    /// One parameter
    bool FitLocationMLE(const std::vector<double> &Sample);
    bool FitScaleMLE(const std::vector<double> &Sample);
    bool FitAsymmetryMLE(const std::vector<double> &Sample);
    /// Two parameters
    bool FitLocationAndScaleMLE(const std::vector<double> &Sample);
    bool FitLocationAndAsymmetryMLE(const std::vector<double> &Sample);
    bool FitScaleAndAsymmetryMLE(const std::vector<double> &Sample);
    /// All parameters
    bool FitMLE(const std::vector<double> &Sample);

private:
    /**
     * @brief GetAsymmetryFromSkewness
     * @param sample
     * @return numeric solution of equation Skewness() = skewness
     */
    double GetAsymmetryFromSkewness(double skewness);
public:
    /// Method of moments
    /// One parameter
    bool FitLocationMM(const std::vector<double> &Sample);
    bool FitScaleMM(const std::vector<double> &Sample);
    bool FitAsymmetryMM(const std::vector<double> &Sample);
    /// Two parameters
    bool FitLocationAndScaleMM(const std::vector<double> &Sample);
    bool FitLocationAndAsymmetryMM(const std::vector<double> &Sample);
    bool FitScaleAndAsymmetryMM(const std::vector<double> &Sample);
    /// All parameters
    bool FitMM(const std::vector<double> &Sample);
};

#endif // LAPLACERAND_H
