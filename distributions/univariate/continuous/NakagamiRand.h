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
    double sigma; /// w / m

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

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
};


/**
 * @brief The ChiRand class
 * Chi distribution
 *
 * Notation: X ~ Chi(n, σ)
 *
 * Related distributions:
 * X ~ Nakagami(n/2, nσ^2)
 *
 * @todo figure out the difference between chi-sigma and nakagami-sigma
 */
class RANDLIBSHARED_EXPORT ChiRand : public NakagamiRand
{
protected:
    double sigma;
    double sigmaSqInv; /// 1.0 / sigma^2

public:
    explicit ChiRand(int degree, double scale = 1.0);
    std::string Name() const override;

private:
    using NakagamiRand::SetParameters;
    using NakagamiRand::GetShape;
    using NakagamiRand::GetSpread;

public:
    void SetParameters(int degree, double scale = 1.0);
    inline int GetDegree() const { return 2 * NakagamiRand::GetShape(); }
    inline double GetScale() const { return sigma; }

private:
    double skewnessImpl(double mean, double sigma) const;

public:
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
 * X ~ Chi(3, σ)
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public ChiRand
{
public:
    explicit MaxwellBoltzmannRand(double scale);
    std::string Name() const override;

    double f(double x) const override;
    double F(double x) const override;

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
 * X ~ Chi(2, σ)
 *
 */
class RANDLIBSHARED_EXPORT RayleighRand : public ChiRand
{
public:
    explicit RayleighRand(double scale = 1);
    std::string Name() const override;

    void SetScale(double scale);

    double f(double x) const override;
    double F(double x) const override;

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
    bool FitScaleMLE(const std::vector<double> &Sample);
    bool FitScaleUMVU(const std::vector<double> &Sample);
};


#endif // NAKAGAMIRAND_H
