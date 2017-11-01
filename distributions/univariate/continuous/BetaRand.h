#ifndef BETARAND_H
#define BETARAND_H

#include "GammaRand.h"
#include "NormalRand.h"

/**
 * @brief The BetaDistribution class <BR>
 * Abstract class for Beta distribution
 *
 * f(x | α, β) = y^{α-1} (1-y)^{β-1} / B(α, β), <BR>
 * where y = (x-a)/(b-a) and B(α, β) denotes Beta function
 *
 * Notation: X ~ B(α, β, a, b) or X ~ B(α, β) for a=0 and b=1
 *
 * Related distributions (a=0, b=1): <BR>
 * 1 − X ~ B(β, α) <BR>
 * X / (1 - X) ~ B'(α, β) <BR>
 * X = Y / (Y + Z), where Y ~ Γ(α) and Z ~ Γ(β) <BR>
 * βX / α(1 - X) ~ F(2α, 2β)
 */
class RANDLIBSHARED_EXPORT BetaDistribution : public ContinuousDistribution
{
protected:
    double alpha = 1; ///< first shape α
    double beta = 1; ///< second shape β
    double a = 0; ///< min bound
    double b = 1; ///< max bound
    double bma = 1; ///< b-a
    double bmaInv = 1; ///< 1/(b-a)
    double logBma = 0; ///< log(b-a)

    GammaRand GammaRV1{}, GammaRV2{};

private:
    static constexpr double edgeForGenerators = 8.0;
    double logBetaFun = 0; ///< log(B(α, β)
    double betaFun = 1; ///< B(α, β)

    /// coefficients for generators
    struct genCoef_t {
        double s, t, u;
    } genCoef = {0, 0, 0};

protected:
    BetaDistribution(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1);
    virtual ~BetaDistribution() {}

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

private:
    enum GENERATOR_ID {
        UNIFORM, ///< standard uniform variate
        ARCSINE, ///< arcsine method
        CHENG, ///< Cheng's method
        REJECTION_UNIFORM, ///< rejection method from uniform distribution for specific value of shapes α = β = 1.5
        REJECTION_UNIFORM_EXTENDED, ///< rejection method from uniform distribution accelerated by using exponential distribution
        REJECTION_NORMAL, ///< rejection method normal distribution
        JOHNK, ///< Johnk's method
        ATKINSON_WHITTAKER, ///< Atkinson-Whittaker's method
        GAMMA_RATIO ///< ratio of two gamma variables
    };

    /**
     * @fn getIdOfUsedGenerator
     * @return id of used variate generator according to the shapes
     */
    GENERATOR_ID getIdOfUsedGenerator() const;

    /**
     * @fn setCoefficientsForGenerator
     */
    void setCoefficientsForGenerator();

protected:
    /**
     * @fn SetShapes
     * @param shape1 α
     * @param shape2 β
     */
    void SetShapes(double shape1, double shape2);
    /**
     * @fn SetSupport
     * @param minValue a
     * @param maxValue b
     */
    void SetSupport(double minValue, double maxValue);

public:
    /**
     * @fn GetAlpha
     * @return first shape α
     */
    inline double GetAlpha() const { return alpha; }
    /**
     * @fn GetBeta
     * @return second shape β
     */
    inline double GetBeta() const { return beta; }
    /**
     * @fn GetBetaFunction
     * @return B(α, β)
     */
    inline double GetBetaFunction() const { return betaFun; }
    /**
     * @fn GetLogBetaFunction
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
     * @fn variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for α = β = 1.5
     */
    double variateRejectionUniform() const;

    /**
     * @fn variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for 1 < α = β < 2 and α != 1.5
     */
    double variateRejectionUniformExtended() const;

    /**
     * @fn variateArcsine
     * Arcsine beta generator
     * @return beta variate for α = β = 0.5
     */
    double variateArcsine() const;

    /**
     * @fn variateRejectionNormal
     * Symmetric beta generator via rejection from the normal density
     * @return beta variate for equal shape parameters > 2
     */
    double variateRejectionNormal() const;

    /**
     * @fn variateJohnk
     * Johnk's beta generator
     * @return beta variate for small shape parameters < 1
     */
    double variateJohnk() const;

    /**
     * @fn variateCheng
     * Cheng's beta generator
     * @return beta variate for max(α, β) > 1 and min(α, β) > 0.5
     */
    double variateCheng() const;

    /**
     * @fn variateAtkinsonWhittaker
     * Atkinson-Whittaker beta generator
     * @return beta variate for max(α, β) < 1 and α + β > 1
     */
    double variateAtkinsonWhittaker() const;

    /**
     * @fn variateGammaRatio
     * Gamma ratio beta generator
     * @return beta variate for the rest variations of shapes
     */
    double variateGammaRatio() const;

public:
    double Mean() const override;
    /**
     * @fn GeometricMean
     * @return E[ln(X)]
     */
    double GeometricMean() const;
    double Variance() const override;
    /**
     * @fn GeometricVariance
     * @return Var(ln(X))
     */
    double GeometricVariance() const;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
};

/**
 * @brief The BetaRand class <BR>
 * Beta distribution
 */
class RANDLIBSHARED_EXPORT BetaRand : public BetaDistribution
{
public:
    BetaRand(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1) : BetaDistribution(shape1, shape2, minValue, maxValue) {}
    String Name() const override;

    using BetaDistribution::SetShapes;
    using BetaDistribution::SetSupport;

    /**
     * @fn FitAlphaMM
     * set α, estimated via method of moments
     * @param sample
     */
    void FitAlphaMM(const std::vector<double> &sample);
    /**
     * @fn FitBetaMM
     * set β, estimated via method of moments
     * @param sample
     */
    void FitBetaMM(const std::vector<double> &sample);
};


/**
 * @brief The ArcsineRand class <BR>
 * Arcsine distribution
 *
 * Notation: X ~ Arcsine(α)
 * 
 * Related distributions: <BR>
 * X ~ B(1 - α, α, 0, 1)
 */
class RANDLIBSHARED_EXPORT ArcsineRand : public BetaDistribution
{
public:
    ArcsineRand(double shape = 0.5, double minValue = 0, double maxValue = 1) : BetaDistribution(1.0 - shape, shape, minValue, maxValue) {}
    String Name() const override;
    void SetShape(double shape);
    inline double GetShape() const { return beta; }
};


/**
 * @brief The BaldingNicholsRand class <BR>
 * Balding-Nichols distribution
 *
 * Notation: X ~ Balding-Nichols(F, p)
 *
 * Related distributions: <BR>
 * X ~ B(p * F', (1 - p) * F', 0, 1) for F' = (1 - F) / F
 */
class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaDistribution
{
    double F = 0.5, p = 0.5;
public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    String Name() const override;

    void SetFixatingIndexAndFrequency(double fixatingIndex, double frequency);
    inline double GetFixatingIndex() const { return F; }
    inline double GetFrequency() const { return p; }
};

#endif // BETARAND_H
