#ifndef BETARAND_H
#define BETARAND_H

#include "ContinuousRand.h"
#include "GammaRand.h"
#include "NormalRand.h"

/**
 * @brief The BetaRand class
 */
class RANDLIBSHARED_EXPORT BetaRand : public ContinuousRand
{

protected:
    GammaRand X, Y;
    NormalRand N; // TODO: we need to storage N OR (X AND Y)
    static constexpr double edgeForGenerators = 8.0;
    double pdfCoef;
    double variateCoef;

public:
    BetaRand(double shape1 = 1, double shape2 = 1);
    virtual std::string name() override;

    void setParameters(double shape1, double shape2);
    void setAlpha(double shape1);
    void setBeta(double shape2);
    inline double getAlpha() const { return X.getShape(); }
    inline double getBeta() const { return Y.getShape(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData);

private:
    double variateForSmallEqualParameters() const;
    double variateForLargeEqualParameters() const;
    double variateForDifferentParameters() const;

    void setVariateConstants();

public:
    double Mean() const override;
    double Variance() const override;
    
    double Quantile(double p) const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    /**
     * @brief getInvBetaFunction
     * @return 1 / B(alpha, beta)
     */
    inline double getInverseBetaFunction() const { return pdfCoef; }
};

#endif // BETARAND_H
