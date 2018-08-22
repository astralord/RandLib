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
template < typename RealType = double >
class RANDLIBSHARED_EXPORT BetaDistribution : public ContinuousDistribution<RealType>
{
protected:
    double alpha = 1; ///< first shape α
    double beta = 1; ///< second shape β
    double a = 0; ///< min bound
    double b = 1; ///< max bound
    double bma = 1; ///< b-a
    double bmaInv = 1; ///< 1/(b-a)
    double logbma = 0; ///< log(b-a)

    GammaRand<RealType> GammaRV1{}, GammaRV2{};

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
    RealType MinValue() const override { return a; }
    RealType MaxValue() const override { return b; }

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
    GENERATOR_ID getIdOfUsedGenerator() const {
        if (alpha < 1 && beta < 1 && alpha + beta > 1)
            return ATKINSON_WHITTAKER;

        if (RandMath::areClose(alpha, beta)) {
            if (RandMath::areClose(alpha, 1.0))
                return UNIFORM;
            else if (RandMath::areClose(alpha, 0.5))
                return ARCSINE;
            else if (RandMath::areClose(alpha, 1.5))
                return REJECTION_UNIFORM;
            else if (alpha > 1)
                return (alpha < 2) ? REJECTION_UNIFORM_EXTENDED : REJECTION_NORMAL;
        }
        if (std::min(alpha, beta) > 0.5 && std::max(alpha, beta) > 1)
            return CHENG;
        return (alpha + beta < 2) ? JOHNK : GAMMA_RATIO;
    }

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

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

private:
    /**
     * @fn variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for α = β = 1.5
     */
    RealType variateRejectionUniform() const;

    /**
     * @fn variateRejectionUniform
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for 1 < α = β < 2 and α != 1.5
     */
    RealType variateRejectionUniformExtended() const;

    /**
     * @fn variateArcsine
     * Arcsine beta generator
     * @return beta variate for α = β = 0.5
     */
    RealType variateArcsine() const;

    /**
     * @fn variateRejectionNormal
     * Symmetric beta generator via rejection from the normal density
     * @return beta variate for equal shape parameters > 2
     */
    RealType variateRejectionNormal() const;

    /**
     * @fn variateJohnk
     * Johnk's beta generator
     * @return beta variate for small shape parameters < 1
     */
    RealType variateJohnk() const;

    /**
     * @fn variateCheng
     * Cheng's beta generator
     * @return beta variate for max(α, β) > 1 and min(α, β) > 0.5
     */
    RealType variateCheng() const;

    /**
     * @fn variateAtkinsonWhittaker
     * Atkinson-Whittaker beta generator
     * @return beta variate for max(α, β) < 1 and α + β > 1
     */
    RealType variateAtkinsonWhittaker() const;

    /**
     * @fn variateGammaRatio
     * Gamma ratio beta generator
     * @return beta variate for the rest variations of shapes
     */
    RealType variateGammaRatio() const;

public:
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    /**
     * @fn GeometricMean
     * @return E[ln(X)]
     */
    long double GeometricMean() const;
    long double Variance() const override;
    /**
     * @fn GeometricVariance
     * @return Var(ln(X))
     */
    long double GeometricVariance() const;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
    /**
     * @brief MeanAbsoluteDeviation
     * @return E[|X - E[X]|]
     */
    long double MeanAbsoluteDeviation() const;

protected:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

    static constexpr char ALPHA_ZERO[] = "Possibly one or more elements of the sample coincide with the lower boundary a.";
    static constexpr char BETA_ZERO[] = "Possibly one or more elements of the sample coincide with the upper boundary b.";
};

/**
 * @brief The BetaRand class <BR>
 * Beta distribution
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT BetaRand : public BetaDistribution<RealType>
{
public:
    BetaRand(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1) : BetaDistribution<RealType>(shape1, shape2, minValue, maxValue) {}
    String Name() const override;

    using BetaDistribution<RealType>::SetShapes;
    using BetaDistribution<RealType>::SetSupport;

public:
    /**
     * @brief GetSampleLogMeanNorm
     * @param sample
     * @return mean average of log(x) for x from normalized sample
     */
    long double GetSampleLogMeanNorm(const std::vector<RealType> &sample) const;
    /**
     * @brief GetSampleLog1pMeanNorm
     * @param sample
     * @return mean average of log(1+x) for x from normalized sample
     */
    long double GetSampleLog1pMeanNorm(const std::vector<RealType> &sample) const;
    /**
     * @brief GetSampleLog1mMeanNorm
     * @param sample
     * @return mean average of log(1-x) for x from normalized sample
     */
    long double GetSampleLog1mMeanNorm(const std::vector<RealType> &sample) const;

    /**
     * @fn FitAlpha
     * set α, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG normalized sample average of ln(X)
     * @param meanNorm normalized sample average
     */
    void FitAlpha(long double lnG, long double meanNorm);

    /**
     * @fn FitAlpha
     * set α, estimated via maximum likelihood
     * @param sample
     */
    void FitAlpha(const std::vector<RealType> &sample);

    /**
     * @fn FitBeta
     * set β, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG1m normalized sample average of ln(1-X)
     * @param meanNorm normalized sample average
     */
    void FitBeta(long double lnG1m, long double meanNorm);

    /**
     * @fn FitBeta
     * set β, estimated via maximum likelihood
     * @param sample
     */
    void FitBeta(const std::vector<RealType> &sample);

    /**
     * @fn FitShapes
     * set α and β, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG sample average of ln(X)
     * @param lnG1m sample average of ln(1-X)
     * @param mean sample average
     * @param variance sample variance
     */
    void FitShapes(long double lnG, long double lnG1m, long double mean, long double variance);

    /**
     * @fn FitShapes
     * set α and β, estimated via maximum likelihood
     * @param sample
     */
    void FitShapes(const std::vector<RealType> &sample);
};


/**
 * @brief The ArcsineRand class <BR>
 * Arcsine distribution
 *
 * Notation: X ~ Arcsine(α, a, b)
 * 
 * Related distributions: <BR>
 * X ~ B(1 - α, α, a, b)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT ArcsineRand : public BetaDistribution<RealType>
{
public:
    ArcsineRand(double shape = 0.5, double minValue = 0, double maxValue = 1) : BetaDistribution<RealType>(1.0 - shape, shape, minValue, maxValue) {}
    String Name() const override;

    using BetaDistribution<RealType>::SetSupport;

    void SetShape(double shape);
    inline double GetShape() const { return this->beta; }

    /**
     * @fn FitShape
     * set α and β, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG average of all ln(X)
     * @param lnG1m average of all ln(1-X)
     */
    void FitShape(long double lnG, long double lnG1m);

    /**
     * @fn FitShape
     * set α, estimated via maximum likelihood
     * @param sample
     */
    void FitShape(const std::vector<RealType> &sample);
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
template < typename RealType = double >
class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaDistribution<RealType>
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
