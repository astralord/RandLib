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
    double alpha, beta;
    GammaRand X, Y;
    double gammaA, gammaB, pdfCoef;
public:
    BetaRand(double shape1, double shape2);
    virtual void setName() override;

    void setParameters(double shape1, double shape2);
    void setAlpha(double shape1);
    void setBeta(double shape2);
    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData);

private:
    double variateForEqualParameters() const;
    double variateForDifferentParameters() const;

public:
    double E() const override { return alpha / (alpha + beta); }
    double Var() const override {
        double denominator = alpha + beta;
        denominator *= denominator * (denominator + 1);
        return alpha * beta / denominator;
    }
};

#endif // BETARAND_H
