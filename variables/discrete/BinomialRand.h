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
    virtual std::string name() override;

    void setNumber(int number);
    inline double getNumber() const { return n; }

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    std::complex<double> CF(double t) const override;

    double E() const override { return n * p; }
    double Var() const override { return n * p * (1 - p); }
};


#endif // BINOMIALRAND_H
