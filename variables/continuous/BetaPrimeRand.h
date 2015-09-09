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

    double E() const override;
    double Var() const override;
    
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // BETAPRIMERAND_H
