#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

#include "GammaRand.h"

/**
 * @brief The NakagamiDistribution class <BR>
 * Abstract class for Nakagami distribution
 *
 * Notation: X ~ Nakagami(m, w)
 *
 * Related distributions: <BR>
 * σX ~ Nakagami(m, wσ^2) <BR>
 * X^2 ~ Γ(m, m / w)
 */
class RANDLIBSHARED_EXPORT NakagamiDistribution : public ContinuousDistribution
{
    double m = 0.5; ///< shape
    double w = 1; ///< spread
    GammaRand Y {};
    double lgammaShapeRatio = 0; ///< log(Γ(m + 0.5) / Γ(m))

protected:
    NakagamiDistribution(double shape = 0.5, double spread = 1);

public:
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

protected:
    /**
     * @fn SetParameters
     * @param shape m
     * @param spread w
     */
    void SetParameters(double shape, double spread);

public:
    /**
     * @fn GetShape
     * @return shape m
     */
    inline double GetShape() const { return m; }
    /**
     * @fn GetSpread
     * @return spread w
     */
    inline double GetSpread() const { return w; }
    /**
     * @fn GetLogGammaFunction
     * @return log(Γ(m))
     */
    inline double GetLogGammaFunction() const { return Y.GetLogGammaShape(); }
    /**
     * @fn GetLogGammaShapeRatio
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
    double Skewness() const override;
    double FourthMoment() const override;
    double ExcessKurtosis() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};


/**
 * @brief The NakagamiRand class <BR>
 * Nakagami distribution
 */
class RANDLIBSHARED_EXPORT NakagamiRand : public NakagamiDistribution
{
public:
    NakagamiRand(double shape = 0.5, double spread = 1) : NakagamiDistribution(shape, spread) {}
    String Name() const override;
    using NakagamiDistribution::SetParameters;
};


/**
 * @brief The ChiRand class <BR>
 * Chi distribution
 *
 * Notation: X ~ χ(k)
 *
 * Related distributions: <BR>
 * X ~ Nakagami(k/2, k) <BR>
 * X^2 ~ χ^2(k) <BR>
 * X^2 ~ Γ(k/2, 0.5)
 */
class RANDLIBSHARED_EXPORT ChiRand : public NakagamiDistribution
{

public:
    explicit ChiRand(int degree);
    String Name() const override;

public:
    /**
     * @fn SetDegree
     * set degree k
     * @param degree
     */
    void SetDegree(int degree);
    /**
     * @fn GetDegree
     * @return degree k
     */
    inline int GetDegree() const { return 2 * NakagamiDistribution::GetShape(); }

    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The MaxwellBoltzmannRand class <BR>
 * Maxwell-Boltzmann distribution
 *
 * Notation: X ~ MB(σ)
 *
 * Related distributions: <BR>
 * X / σ ~ χ(3) <BR>
 * X ~ Nakagami(1.5, 3σ^2)
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public NakagamiDistribution
{
    double sigma = 1; ///< scale σ
public:
    explicit MaxwellBoltzmannRand(double scale);
    String Name() const override;

public:
    /**
     * @fn SetScale
     * set scale σ
     * @param scale
     */
    void SetScale(double scale);
    /**
     * @fn GetScale
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
 * @brief The RayleighRand class <BR>
 * Rayleigh distribution
 *
 * Notation: X ~ Rayleigh(σ)
 *
 * Related distributions:
 * X / σ ~ χ(2)
 * X ~ Nakagami(1, 2σ^2)
 */
class RANDLIBSHARED_EXPORT RayleighRand : public NakagamiDistribution
{
    double sigma = 1; ///< scale σ
public:
    explicit RayleighRand(double scale = 1);
    String Name() const override;

public:
    /**
     * @fn SetScale
     * set scale σ
     * @param scale
     */
    void SetScale(double scale);
    /**
     * @fn GetScale
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
     * @fn FitScale
     * set scale via maximum-likelihood method if unbiased = false,
     * otherwise set scale, returned by uniformly minimum variance unbiased estimator
     * @param sample
     */
    void FitScale(const std::vector<double> &sample, bool unbiased = false);
};


#endif // NAKAGAMIRAND_H
