#ifndef NEGATIVEBINOMIALRAND_H
#define NEGATIVEBINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"
#include "../continuous/GammaRand.h"


/**
 * @brief The NegativeBinomialRand class
 * Negative binomial distribution
 * X ~ NB(r, p)
 *
 * P(X = k) = binom(k + r - 1, k) p^r (1-p)^k
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
    double logQ;
    static constexpr int tableSize = 16;
    double table[tableSize];

public:
    NegativeBinomialRand(T number, double probability);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

private:
    void setValidParameters(T number, double probability);
public:
    void setParameters(T number, double probability);
    inline double getProbability() const { return p; }
    inline T getNumber() const { return r; }

    double P(int k) const override;
    double F(int k) const override;
protected:

    enum GENERATOR_ID {
        TABLE,
        EXPONENTIAL,
        GAMMA_POISSON
    };

    GENERATOR_ID getIdOfUsedGenerator() const;

    int variateGeometricByTable() const;
    int variateGeometricThroughExponential() const;
private:
    int variateByTable() const;
    int variateThroughExponential() const;
    int variateThroughGammaPoisson() const;

public:
    int variate() const override;
    void sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    
    std::complex<double> CF(double t) const override;
    
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


typedef NegativeBinomialRand<int> PascalRand;
typedef NegativeBinomialRand<double> PolyaRand;

#endif // NEGATIVEBINOMIALRAND_H
