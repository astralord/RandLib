#ifndef NEGATIVEBINOMIALRAND_H
#define NEGATIVEBINOMIALRAND_H

#include "DiscreteRand.h"
#include "GeometricRand.h"
#include "PoissonRand.h"
#include "../continuous/GammaRand.h"

template < typename T >
class RANDLIBSHARED_EXPORT NegativeBinomialRand : public DiscreteRand
{
    double p;
    T r;

    double pdfCoef;

    GeometricRand G;
    GammaRand Y;

public:
    NegativeBinomialRand(T number, double probability);
    virtual std::string name() override;

    void setParameters(T number, double probability);
    inline double getProbability() const { return p; }
    inline T getNumber() const { return r; }

    double P(int r) const override;
    double F(double x) const override;
    double variate() const override;

private:
    double variateThroughGeometric() const;
    double variateThroughGammaPoisson() const;

public:
    double Mean() const override;
    double Variance() const override;
    
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

typedef NegativeBinomialRand<int> PascaleRand;
typedef NegativeBinomialRand<double> PolyaRand;


#endif // NEGATIVEBINOMIALRAND_H
