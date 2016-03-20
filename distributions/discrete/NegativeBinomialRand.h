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

private:
    T r;
    GammaRand Y;
    double pdfCoef;
    static constexpr int tableSize = 16;
    double table[tableSize];

public:
    NegativeBinomialRand(T number, double probability);
    std::string name() override;

private:
    void setValidParameters(T number, double probability);
public:
    void setParameters(T number, double probability);
    inline double getProbability() const { return p; }
    inline T getNumber() const { return r; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(QVector<double> &outputData) const override;

protected:
    double variateGeometricByTable() const;
    double variateGeometricThroughExponential() const;
private:
    double variateThroughGeometric() const;
    double variateThroughGammaPoisson() const;

public:
    double Mean() const override;
    double Variance() const override;
    
    std::complex<double> CF(double t) const override;
    
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


typedef NegativeBinomialRand<int> PascalRand;
typedef NegativeBinomialRand<double> PolyaRand;

#endif // NEGATIVEBINOMIALRAND_H
