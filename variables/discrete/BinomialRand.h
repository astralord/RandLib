#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "BernoulliRand.h"

/**
 * @brief The BinomialRand class
 */
class RANDLIBSHARED_EXPORT BinomialRand : public DiscreteRand
{
    int n;
    double p;

    BernoulliRand B;

public:
    BinomialRand(int number, double probability);
    std::string name() override;

    void setNumber(int number);
    inline double getNumber() const { return n; }

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override { return n * p; }
    double Variance() const override { return n * p * (1 - p); }

    std::complex<double> CF(double t) const override;
    
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


#endif // BINOMIALRAND_H
