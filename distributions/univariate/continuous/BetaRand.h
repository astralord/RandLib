#ifndef BETARAND_H
#define BETARAND_H

#include "GammaRand.h"
#include "NormalRand.h"

/**
 * @brief The BetaRand class
 * Beta distribution
 * X ~ Beta(\alpha, \beta)
 *
 * f(x|\alpha, \beta) = x^{\alpha-1} (1-x)^{\beta-1} / B(\alpha, \beta)
 */
class RANDLIBSHARED_EXPORT BetaRand : public ContinuousDistribution
{
protected:
    double alpha, beta; /// hashed parameters
    double a, b, bma;

private:
    GammaRand X, Y;
    NormalRand N; // TODO: we need to storage N OR (X AND Y)
    static constexpr double edgeForGenerators = 8.0;
    double pdfCoef, cdfCoef;
    double variateCoef;

public:
    BetaRand(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1);
    virtual ~BetaRand() {}
    std::string name() override;

    void setShapes(double shape1, double shape2);
    void setSupport(double minValue, double maxValue);
    void setAlpha(double shape1);
    void setBeta(double shape2);
    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }
    inline double getMin() const { return a; }
    inline double getMax() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData) const override;

private:
    double variateArcsine() const;                  /// alpha = beta = 0.5
    double variateForSmallEqualParameters() const;  /// alpha = beta, 1 < alpha < edgeForGenerators
    double variateForLargeEqualParameters() const;  /// alpha = beta, alpha > edgeForGenerators
    double variateForDifferentParameters() const;   /// otherwise

    void setVariateConstants();

public:
    double Mean() const override;
    double Variance() const override;
    
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    /**
     * @brief getInvBetaFunction
     * @return 1 / B(alpha, beta)
     */
    inline double getInverseBetaFunction() const { return cdfCoef; }
};


/**
 * @brief The ArcsineRand class
 */
class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
public:
    ArcsineRand(double shape = 0.5, double minValue = 0, double maxValue = 1);
    std::string name() override;

    void setShape(double shape);
    inline double getShape() const { return beta; }

protected:
    /// prohibit to use beta's getters and setters
    using BetaRand::setShapes;
    using BetaRand::setAlpha;
    using BetaRand::setBeta;
};


/**
 * @brief The BaldingNicholsRand class
 */
class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaRand
{
    double p, F;
public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    std::string name() override;

    void setParameters(double fixatingIndex, double frequency);
    inline double getFrequency() { return p; }
    inline double getFixatingIndex() { return F; }

private:
    using BetaRand::setShapes;
    using BetaRand::setAlpha;
    using BetaRand::setBeta;
};

#endif // BETARAND_H
