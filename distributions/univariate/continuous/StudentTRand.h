#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "NakagamiRand.h"

/**
 * @brief The StudentTRand class <BR>
 * Student's t-distribution
 *
 * Notation: X ~ t(ν, μ, σ)
 * If X ~ t(1, μ, σ), then X ~ Cauchy(μ, σ)
 * X -> Normal(μ, σ) for t -> ∞
 */
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousDistribution
{
    double nu = 1; ///< degree ν
    double mu = 0; ///< location μ
    double sigma = 1; ///< scale σ
    double logSigma = 0; ///< log(σ)
    NakagamiRand Y{};
    double pdfCoef = -M_LNPI; ///< coefficient for faster pdf calculation
    double nup1Half = 1; ///< 0.5 * (ν + 1)
    double logBetaFun = M_LNPI; ///< log(B(0.5 * ν, 0.5))

public:
    explicit StudentTRand(double degree = 1.0, double location = 0.0, double scale = 1.0);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetDegree(double degree);
    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetDegree() const { return nu; }
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return sigma; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // STUDENTTRAND_H
