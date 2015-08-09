#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "DiscreteRand.h"
#include "BernoulliRand.h"

/**
 * @brief The BinomialRand class
 */
class RANDLIBSHARED_EXPORT BinomialRand : public DiscreteRand<int>
{
    int n;
    double p;

    BernoulliRand B;

public:
    BinomialRand(int number, double probability);

    void setN(int number);
    inline double getN() const { return n; }

    void setP(double probability);
    inline double getP() const { return p; }

    virtual double P(int k) const override;
    virtual double F(double x) const override;
    virtual double variate() const override;

    // WRITE var!
    double E() const override { return n * p; }
    double Var() const override { return 1; }
};


#endif // BINOMIALRAND_H
