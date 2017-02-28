#ifndef NEGATIVEBINOMIALRAND_H
#define NEGATIVEBINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"
#include "../continuous/GammaRand.h"


/**
 * @brief The NegativeBinomialRand class
 * Negative binomial distribution
 *
 * P(X = k) = C(k + r - 1, k) p^r (1-p)^k
 *
 * Notation: X ~ NB(r, p)
 *
 * Related distributions:
 * If X ~ NB(1, p), then X ~ Geometric(p)
 */
template < typename T >
class RANDLIBSHARED_EXPORT NegativeBinomialRand : public DiscreteDistribution
{
protected:
    double p, q;
    double logProb, log1mProb; /// log(p) and log(q)

private:
    T r;
    double pdfCoef;
    double qDivP; /// q / p
    static constexpr int tableSize = 16;
    double table[tableSize];
    GammaRand GammaRV;

public:
    NegativeBinomialRand(T number, double probability);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return INT_MAX; }

private:
    void SetValidParameters(T number, double probability);
public:
    void SetParameters(T number, double probability);
    inline double GetProbability() const { return p; }
    inline T GetNumber() const { return r; }

    double P(int k) const override;
    double logP(int k) const override;
    double F(int k) const override;
    double S(int k) const override;

protected:
    enum GENERATOR_ID {
        TABLE,
        EXPONENTIAL,
        GAMMA_POISSON
    };

    GENERATOR_ID GetIdOfUsedGenerator() const;

    int variateGeometricByTable() const;
    int variateGeometricThroughExponential() const;
private:
    int variateByTable() const;
    int variateThroughExponential() const;
    int variateThroughGammaPoisson() const;

public:
    int Variate() const override;
    void Sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /// Method of moments
    bool FitNumberMM(const std::vector<int> &sample);
    bool FitProbabilityMM(const std::vector<int> & sample);
    bool FitMM(const std::vector<int> & sample);

    /// Maximum likelihood
    bool FitMLE(const std::vector<int> &sample);
};


typedef NegativeBinomialRand<int> PascalRand;
typedef NegativeBinomialRand<double> PolyaRand;

#endif // NEGATIVEBINOMIALRAND_H
