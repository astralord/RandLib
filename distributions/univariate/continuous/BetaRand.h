#ifndef BETARAND_H
#define BETARAND_H

#include "GammaRand.h"
#include "NormalRand.h"

/**
 * @brief The BetaRand class
 * Beta distribution
 *
 * f(x | α, β) = x^{α-1} (1-x)^{β-1} / B(α, β),
 * where B(α, β) denotes Beta function
 *
 * Notation: X ~ Beta(α, β)
 *
 * Related distributions:
 * 1 − X ~ Beta(β, α)
 * X / (1 - X) ~ Beta-Prime(α, β)
 * X = Y / (Y + Z), where Y ~ Gamma(α) and Z ~ Gamma(β)
 * βX / α(1 - X) ~ F(2α, 2β)
 */
class RANDLIBSHARED_EXPORT BetaRand : public ContinuousDistribution
{
protected:
    /// parameters of distribution
    double alpha, beta;
    double a, b, bma, bmaInv, logBma;

    GammaRand GammaRV1, GammaRV2;

private:
    static constexpr double edgeForGenerators = 8.0;
    double logBetaFun; /// log(B(α, β)
    double betaFun; /// B(α, β)

    /// coefficients for generators
    double s, t, u;

public:
    BetaRand(double shape1 = 1, double shape2 = 1);
    BetaRand(double shape1, double shape2, double minValue, double maxValue);
    virtual ~BetaRand() {}
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

private:

    enum GENERATOR_ID {
        UNIFORM,
        ARCSINE,
        CHENG,
        REJECTION_UNIFORM,
        REJECTION_UNIFORM_EXTENDED,
        REJECTION_NORMAL,
        JOHNK,
        ATKINSON_WHITTAKER,
        GAMMA_RATIO
    };

    /**
     * @brief getIdOfUsedGenerator
     * @return id of used variate generator according to the shapes
     */
    GENERATOR_ID getIdOfUsedGenerator() const;

    /**
     * @brief setCoefficientsForGenerator
     */
    void setCoefficientsForGenerator();

public:
    /**
     * @brief SetShapes
     * @param shape1 α
     * @param shape2 β
     */
    void SetShapes(double shape1, double shape2);
    /**
     * @brief SetSupport
     * @param minValue a
     * @param maxValue b
     */
    void SetSupport(double minValue, double maxValue);
    /**
     * @brief GetAlpha
     * @return α
     */
    inline double GetAlpha() const { return alpha; }
    /**
     * @brief GetBeta
     * @return β
     */
    inline double GetBeta() const { return beta; }
    /**
     * @brief GetBetaFunction
     * @return B(α, β)
     */
    inline double GetBetaFunction() const { return betaFun; }
    /**
     * @brief GetLogBetaFunction
     * @return log(B(α, β))
     */
    inline double GetLogBetaFunction() const { return logBetaFun; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

private:
    /**
     * @brief variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for α = β = 1.5
     */
    double variateRejectionUniform() const;

    /**
     * @brief variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for 1 < α = β < 2 and α != 1.5
     */
    double variateRejectionUniformExtended() const;

    /**
     * @brief variateArcsine
     * Arcsine beta generator
     * @return beta variate for α = β = 0.5
     */
    double variateArcsine() const;

    /**
     * @brief variateRejectionNormal
     * Symmetric beta generator via rejection from the normal density
     * @return beta variate for equal shape parameters > 2
     */
    double variateRejectionNormal() const;

    /**
     * @brief variateJohnk
     * Johnk's beta generator
     * @return beta variate for small shape parameters < 1
     */
    double variateJohnk() const;

    /**
     * @brief variateCheng
     * Cheng's beta generator
     * @return beta variate for max(α, β) > 1 and min(α, β) > 0.5
     */
    double variateCheng() const;

    /**
     * @brief variateAtkinsonWhittaker
     * Atkinson-Whittaker beta generator
     * @return beta variate for max(α, β) < 1 and α + β > 1
     */
    double variateAtkinsonWhittaker() const;

    /**
     * @brief variateGammaRatio
     * Gamma ratio beta generator
     * @return beta variate for the rest variations of shapes
     */
    double variateGammaRatio() const;

public:
    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @brief FitAlphaMM
     * set α, estimated via method of moments
     * @param sample
     */
    void FitAlphaMM(const std::vector<double> &sample);
    /**
     * @brief FitBetaMM
     * set β, estimated via method of moments
     * @param sample
     */
    void FitBetaMM(const std::vector<double> &sample);
};


/**
 * @brief The ArcsineRand class
 * Arcsine distribution
 *
 * Notation: X ~ Arcsine(α)
 * 
 * Related distributions
 * X ~ Beta(1 - α, α)
 */
class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
public:
    ArcsineRand(double shape = 0.5, double minValue = 0, double maxValue = 1);
    std::string Name() const override;

    void SetShape(double shape);
    inline double GetShape() const { return beta; }

protected:
    /// prohibit to use beta's setters
    using BetaRand::SetShapes;
    void SetSupport(double minValue, double maxValue) = delete;
    /// and estimators
    void FitAlphaMM(const std::vector<double> &sample) = delete;
    void FitBetaMM(const std::vector<double> &sample) = delete;
};


/**
 * @brief The BaldingNicholsRand class
 * Balding-Nichols distribution
 *
 * Notation: X ~ Balding-Nichols(F, p)
 *
 * Related distributions
 * X ~ Beta(p * F', (1 - p) * F') for F' = (1 - F) / F
 */
class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaRand
{
    double p, F;
public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    std::string Name() const override;

    void SetFixatingIndexAndFrequency(double fixatingIndex, double frequency);
    inline double GetFixatingIndex() const { return F; }
    inline double GetFrequency() const { return p; }

private:
    /// prohibit to use beta's setters
    using BetaRand::SetShapes;
    void SetSupport(double minValue, double maxValue) = delete;
    /// and estimators
    void FitAlphaMM(const std::vector<double> &sample) = delete;
    void FitBetaMM(const std::vector<double> &sample) = delete;
};

#endif // BETARAND_H
