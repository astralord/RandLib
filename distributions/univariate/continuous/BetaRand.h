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
    double a, b, bma;

    GammaRand GammaRV1, GammaRV2;

private:
    static constexpr double edgeForGenerators = 8.0;
    double mLogBetaFun; /// log(Beta(α, β)
    double betaFunInv; /// 1 / Beta(α, β)

    /// coefficients for generators
    double s, t, u;

public:
    BetaRand(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1);
    virtual ~BetaRand() {}
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

    void SetParameters(double shape1, double shape2, double minValue = 0, double maxValue = 1);
    inline double GetAlpha() const { return alpha; }
    inline double GetBeta() const { return beta; }
    inline double GetMin() const { return a; }
    inline double GetMax() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    void Sample(std::vector<double> &outputData) const override;

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
     * @brief GetIdOfUsedGenerator
     * @return id of the used variate generator according to shape parameters
     */
    GENERATOR_ID GetIdOfUsedGenerator() const;

    /**
     * @brief SetCoefficientsForGenerator
     */
    void SetCoefficientsForGenerator();

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

    std::complex<double> CF(double t) const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    /**
     * @brief GetInverseBetaFunction
     * @return 1 / Beta(α, β)
     */
    inline double GetInverseBetaFunction() const { return betaFunInv; }

    /**
     * @brief GetLogBetaFunction
     * @return log Beta(α, β)
     */
    inline double GetLogBetaFunction() const { return -mLogBetaFun; }
};


/**
 * @brief The ArcsineRand class 
 *
 * X ~ Arcsine(α)
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
    /// prohibit to use beta's Getters and Setters
    using BetaRand::SetParameters;
};


/**
 * @brief The BaldingNicholsRand class
 */
class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaRand
{
    double p, F;
public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    std::string Name() const override;

    void SetFixatingIndexAndFrequency(double fixatingIndex, double frequency);
    inline double GetFrequency() const { return p; }
    inline double GetFixatingIndex() const { return F; }

private:
    using BetaRand::SetParameters;
};

#endif // BETARAND_H
