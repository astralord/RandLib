#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

#include "GammaRand.h"

/**
 * @brief The NakagamiRand class
 * Nakagami distribution
 *
 * Notation: X ~ Nakagami(m, w)
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

public:
    double GetLogGammaFunction() const { return Y.GetLogGammaFunction(); }
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
    void SetParameters(int degree);
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
 * X / σ ~ Nakagami(1.5, 3)
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public ChiRand
{
    double sigma;
public:
    explicit MaxwellBoltzmannRand(double scale);
    std::string Name() const override;
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

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};


/**
 * @brief The RayleighRand class
 * Rayleigh distribution
 *
 * Notation: X ~ Rayleigh(σ)
 *
 * Related distributions:
 * X / σ ~ Chi(2)
 * X / σ ~ Nakagami(1, 2)
 */
class RANDLIBSHARED_EXPORT RayleighRand : public ChiRand
{
    double sigma;
public:
    explicit RayleighRand(double scale = 1);
    std::string Name() const override;
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
    void FitScaleMLE(const std::vector<double> &sample);
    void FitScaleUMVU(const std::vector<double> &sample);
};


#endif // NAKAGAMIRAND_H
