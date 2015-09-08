#ifndef BETAPRIMERAND_H
#define BETAPRIMERAND_H

#include "BetaRand.h"

/**
 * @brief The BetaRand class
 */
class RANDLIBSHARED_EXPORT BetaPrimeRand : public BetaRand
{
public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    virtual std::string name() override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData);

    double E() const override { return (Y.getShape() > 1) ? X.getShape() / (Y.getShape() - 1) : INFINITY; }
    double Var() const override {
        double alpha = X.getShape();
        double beta = Y.getShape();
        if (beta <= 2)
            return INFINITY;
        double betaAdj = beta - 1;
        double numerator = alpha * (alpha + betaAdj);
        double denominator = (betaAdj - 1) * betaAdj * betaAdj;
        return numerator / denominator;
    }
    
    inline double Mode() { return (X.getShape() < 1) ? 0 : (X.getShape() - 1) / (Y.getShape() + 1);
    double Skewness() {
        double alpha = X.getShape();
        double beta = Y.getShape();
        if (beta <= 3)
            return INFINITY;
        double aux = alpha + beta - 1;
        double skewness = (beta - 2) / (alpha * aux);
        skewness = std::sqrt(skewness);
        aux += alpha;
        aux += aux;
        return aux * skewness / (beta - 3);
    }
};

#endif // BETAPRIMERAND_H
