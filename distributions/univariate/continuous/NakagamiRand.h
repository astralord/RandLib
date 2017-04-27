#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

#include "GammaRand.h"

/**
 * @brief The NakagamiDistribution class
 * Abstract class for Nakagami distribution
 *
 * Notation: X ~ Nakagami(m, w)
 *
 * Related distributions:
 * σX ~ Nakagami(m, wσ^2)
 * X^2 ~ Gamma(m, m / w)
 */
class RANDLIBSHARED_EXPORT NakagamiDistribution : public ContinuousDistribution
{
    double m, w;
    GammaRand Y;
    /// log(Γ(m + 0.5) / Γ(m))
    double lgammaShapeRatio;

public:
    NakagamiDistribution(double shape = 0.5, double spread = 1);

    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

protected:
    /**
     * @brief SetParameters
     * @param shape m
     * @param spread w
     */
    void SetParameters(double shape, double spread);

public:
    /**
     * @brief GetShape
     * @return shape m
     */
    inline double GetShape() const { return m; }
    /**
     * @brief GetSpread
     * @return spread w
     */
    inline double GetSpread() const { return w; }
    /**
     * @brief GetLogGammaFunction
     * @return log(Γ(m))
     */
    inline double GetLogGammaFunction() const { return Y.GetLogGammaFunction(); }
    /**
     * @brief GetLogGammaShapeRatio
     * @return log(Γ(m + 0.5) / Γ(m))
     */
    inline double GetLogGammaShapeRatio() const { return lgammaShapeRatio; }

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
 * @brief The NakagamiRand class
 * Nakagami distribution
 */
class RANDLIBSHARED_EXPORT NakagamiRand : public NakagamiDistribution
{
public:
    NakagamiRand(double shape = 0.5, double spread = 1) : NakagamiDistribution(shape, spread) {}
    std::string Name() const override;
    using NakagamiDistribution::SetParameters;
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
class RANDLIBSHARED_EXPORT ChiRand : public NakagamiDistribution
{

public:
    explicit ChiRand(int degree);
    std::string Name() const override;

private:
    using NakagamiDistribution::SetParameters;

public:
    void SetDegree(int degree);
    /**
     * @brief GetDegree
     * @return degree k
     */
    inline int GetDegree() const { return 2 * NakagamiDistribution::GetShape(); }

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
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public NakagamiDistribution
{
    double sigma;
public:
    explicit MaxwellBoltzmannRand(double scale);
    std::string Name() const override;

private:
    using NakagamiDistribution::SetParameters;

public:
    void SetScale(double scale);
    /**
     * @brief GetScale
     * @return scale σ
     */
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
class RANDLIBSHARED_EXPORT RayleighRand : public NakagamiDistribution
{
    double sigma;
public:
    explicit RayleighRand(double scale = 1);
    std::string Name() const override;

private:
    using NakagamiDistribution::SetParameters;

public:
    void SetScale(double scale);
    /**
     * @brief GetScale
     * @return scale σ
     */
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
