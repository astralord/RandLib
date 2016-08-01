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
    double alpha, beta; /// hashed parameters
    double a, b, bma;

private:
    static constexpr double edgeForGenerators = 8.0;
    double mLogBeta; /// log(Beta(α, β)
    double betaInv; /// 1 / Beta(α, β)

    /// coefficients for generator:
    double s, t, u;

public:
    BetaRand(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1);
    virtual ~BetaRand() {}
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

    void setParameters(double shape1, double shape2, double minValue = 0, double maxValue = 1);
    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }
    inline double getMin() const { return a; }
    inline double getMax() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(std::vector<double> &outputData) const override;

private:

    enum GENERATOR_ID {
        UNIFORM,
        ARCSINE,
        CHENG,
        REJECTION_UNIFORM,
        REJECTION_NORMAL,
        JOHNK,
        ATKINSON_WHITTAKER
    };

    /**
     * @brief getIdOfUsedGenerator
     * @return id of the used variate generator according to shape parameters
     */
    GENERATOR_ID getIdOfUsedGenerator() const;

    /**
     * @brief setCoefficientsForGenerator
     */
    void setCoefficientsForGenerator();

    /**
     * @brief variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for 1 < α = β < 2
     */
    double variateRejectionUniform() const;

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
     * @return
     */
    double variateCheng() const;

    /**
     * @brief variateAtkinsonWittaker
     * @return
     */
    double variateAtkinsonWhittaker() const;

public:
    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    /**
     * @brief getInvBetaFunction
     * @return 1 / Beta(α, β)
     */
    inline double getInverseBetaFunction() const { return betaInv; }

    /**
     * @brief getLogBetaFunction
     * @return log Beta(α, β)
     */
    inline double getLogBetaFunction() const { return -mLogBeta; }
};


/**
 * @brief The ArcsineRand class
 */
class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
public:
    ArcsineRand(double shape = 0.5, double minValue = 0, double maxValue = 1);
    std::string name() const override;

    void setShape(double shape);
    inline double getShape() const { return beta; }

protected:
    /// prohibit to use beta's getters and setters
    using BetaRand::setParameters;
};


/**
 * @brief The BaldingNicholsRand class
 */
class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaRand
{
    double p, F;
public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    std::string name() const override;

    void setFixatingIndexAndFrequency(double fixatingIndex, double frequency);
    inline double getFrequency() const { return p; }
    inline double getFixatingIndex() const { return F; }

private:
    using BetaRand::setParameters;
};

#endif // BETARAND_H
