#ifndef NEGATIVEBINOMIALRAND_H
#define NEGATIVEBINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"
#include "../continuous/GammaRand.h"


/**
 * @brief The NegativeBinomialRand class
 *
 * P(X = k) = binom(k + r - 1, k) p^r (1-p)^k
 *
 * X ~ NB(r, p)
 */
template < typename T >
class RANDLIBSHARED_EXPORT NegativeBinomialRand : public DiscreteDistribution
{
protected:
    double p, q;

protected:
    T r;
    GammaRand Y;
    double pdfCoef;


public:
    NegativeBinomialRand(T number, double probability);
    std::string name() override;

    virtual void setParameters(T number, double probability);
    inline double getProbability() const { return p; }
    inline T getNumber() const { return r; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

protected:
    double variateThroughGammaPoisson() const;

public:
    double Mean() const override;
    double Variance() const override;
    
    std::complex<double> CF(double t) const override;
    
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


typedef NegativeBinomialRand<int> NegativeBinomialIntRand;
typedef NegativeBinomialRand<double> NegativeBinomialDoubleRand;


/**
 * @brief The PascalRand class
 */
class RANDLIBSHARED_EXPORT PascalRand : public NegativeBinomialIntRand {
    static constexpr int tableSize = 16;
    double table[tableSize];
public:
    PascalRand(int number, double probability);
    std::string name() override;

    void setParameters(int number, double probability) override;

    double P(int k) const override;
    double variate() const override;

private:
    double variateThroughGeometric() const;

protected:
    double variateGeometricByTable() const;
    double variateGeometricThroughExponential() const;
};

/**
 * @brief The PolyaRand class
 */
class RANDLIBSHARED_EXPORT PolyaRand : public NegativeBinomialDoubleRand {
public:
    PolyaRand(double number, double probability);
    std::string name() override;
    void setParameters(double number, double probability) override;
};


#endif // NEGATIVEBINOMIALRAND_H
