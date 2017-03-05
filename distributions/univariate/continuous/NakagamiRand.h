#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

#include "GammaRand.h"

/**
 * @brief The NakagamiRand class
 * Nakagami distribution
 *
 * Notation: X ~ Nakagami(m, w)
 *
 * Related distributions:
 * σX ~ Nakagami(m, wσ^2)
 * X^2 ~ Gamma(m, m / w)
 */
class RANDLIBSHARED_EXPORT NakagamiRand : public ContinuousDistribution
{
    double m, w;
    GammaRand Y;

public:
    NakagamiRand(double shape = 0.5, double spread = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double spread);
    inline double GetShape() const { return m; }
    inline double GetSpread() const { return w; }
    inline double GetLogGammaFunction() const { return Y.GetLogGammaFunction(); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};


/**
 * @brief The ChiRand class
 * Chi distribution
 *
 * Notation: X ~ Chi(k)
 *
 * Related distributions:
 * X ~ Nakagami(k/2, k)
 * X^2 ~ Chi-Squared(k)
 * X^2 ~ Gamma(k/2, 0.5)
 */
class RANDLIBSHARED_EXPORT ChiRand : public NakagamiRand
{

public:
    explicit ChiRand(int degree);
    std::string Name() const override;

private:
    using NakagamiRand::SetParameters;

public:
    void SetDegree(int degree);
    inline int GetDegree() const { return 2 * NakagamiRand::GetShape(); }

    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The MaxwellBoltzmannRand class
 * Maxwell-Boltzmann distribution
 *
 * Notation: X ~ MB(σ)
 *
 * Related distributions:
 * X / σ ~ Chi(3)
 * X ~ Nakagami(1.5, 3σ^2)
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public NakagamiRand
{
    double sigma;
public:
    explicit MaxwellBoltzmannRand(double scale);
    std::string Name() const override;

private:
    using NakagamiRand::SetParameters;

public:
    void SetScale(double scale);
    double GetScale() const { return sigma; }

    double f(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The RayleighRand class
 * Rayleigh distribution
 *
 * Notation: X ~ Rayleigh(σ)
 *
 * Related distributions:
 * X / σ ~ Chi(2)
 * X ~ Nakagami(1, 2σ^2)
 */
class RANDLIBSHARED_EXPORT RayleighRand : public NakagamiRand
{
    double sigma;
public:
    explicit RayleighRand(double scale = 1);
    std::string Name() const override;

private:
    using NakagamiRand::SetParameters;

public:
    void SetScale(double scale);
    double GetScale() const { return sigma; }

    double f(const double & x) const override;
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

public:
    /**
     * @brief FitScaleMLE
     * set scale via maximum-likelihood method
     * @param sample
     */
    void FitScaleMLE(const std::vector<double> &sample);
    /**
     * @brief FitScaleUMVU
     * set scale, returned by uniformly minimum variance unbiased estimator
     * @param sample
     */
    void FitScaleUMVU(const std::vector<double> &sample);
};


#endif // NAKAGAMIRAND_H
