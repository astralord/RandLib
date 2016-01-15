#ifndef BETARAND_H
#define BETARAND_H

#include "GammaRand.h"
#include "NormalRand.h"

/**
 * @brief The BetaRand class
 *
 *
 * X ~ Beta(alpha, beta)
 */
class RANDLIBSHARED_EXPORT BetaRand : public ContinuousRand
{

protected:
    double alpha, beta; /// hashed parameters

private:
    GammaRand X, Y;
    NormalRand N; // TODO: we need to storage N OR (X AND Y)
    static constexpr double edgeForGenerators = 8.0;
    double pdfCoef;
    double variateCoef;

public:
    BetaRand(double shape1 = 1, double shape2 = 1);
    virtual ~BetaRand() {}
    std::string name() override;

    void setParameters(double shape1, double shape2);
    void setAlpha(double shape1);
    void setBeta(double shape2);
    inline double getAlpha() const { return X.getShape(); }
    inline double getBeta() const { return Y.getShape(); }

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
    inline double getInverseBetaFunction() const { return pdfCoef; }
};



class RANDLIBSHARED_EXPORT BaldingNicholsRand : public BetaRand
{
    double p, F;
public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    std::string name() override;

    void setParameters(double fixatingIndex, double frequency);
    inline double getFrequency() { return p; }
    inline double getFixatingIndex() { return F; }

protected:
    using BetaRand::setAlpha;
    using BetaRand::setBeta;
};

#endif // BETARAND_H
