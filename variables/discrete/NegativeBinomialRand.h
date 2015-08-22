#ifndef NEGATIVEBINOMIALRAND_H
#define NEGATIVEBINOMIALRAND_H

#include "DiscreteRand.h"
#include "GeometricRand.h"

template < typename T >
class RANDLIBSHARED_EXPORT NegativeBinomialRand : public DiscreteRand<int>
{
    double r;
    T k;

    GeometricRand G;
public:
    NegativeBinomialRand(T number, double probability);
    virtual void setName() override;

    void setProbability(double probability);
    inline double getProbability() const { return r; }

    void setNumber(T number);
    inline T getNumber() const { return k; }

    //TODO!
    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    //TODO!
    double E() const override { return 0; }
    double Var() const override { return 0; }
};

typedef NegativeBinomialRand<int> PascaleRand;
typedef NegativeBinomialRand<double> PolyaRand;


#endif // NEGATIVEBINOMIALRAND_H
