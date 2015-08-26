#ifndef BETARAND_H
#define BETARAND_H

#include "ContinuousRand.h"
#include "GammaRand.h"

/**
 * @brief The BetaRand class
 */
class RANDLIBSHARED_EXPORT BetaRand : public ContinuousRand
{

protected:
    GammaRand X, Y;
    double pdfCoef;

public:
    BetaRand(double shape1, double shape2);
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
    double variateForEqualParameters() const;
    double variateForDifferentParameters() const;

public:
    double E() const override { return X.getShape() / (X.getShape() + Y.getShape()); }
    double Var() const override {
        double alpha = X.getShape();
        double beta = Y.getShape();
        double denominator = alpha + beta;
        denominator *= denominator * (denominator + 1);
        return alpha * beta / denominator;
    }

    /**
     * @brief getInvBetaFunction
     * @return 1 / B(alpha, beta)
     */
    inline double getInverseBetaFunction() { return pdfCoef; }
};

#endif // BETARAND_H
