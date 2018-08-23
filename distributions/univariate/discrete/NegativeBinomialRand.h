#ifndef NEGATIVEBINOMIALRAND_H
#define NEGATIVEBINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The NegativeBinomialDistribution class <BR>
 * Abstract class for Negative binomial distribution
 *
 * P(X = k) = C(k + r - 1, k) p^r (1-p)^k
 *
 * Notation: X ~ NB(r, p)
 *
 * Related distributions: <BR>
 * If X ~ NB(1, p), then X ~ Geometric(p)
 * If Y ~ Γ(r, p / (1 - p), then Po(Y) ~ NB(r, p)
 */
template< typename IntType = int, typename T = double>
class RANDLIBSHARED_EXPORT NegativeBinomialDistribution : public DiscreteDistribution<IntType>
{
protected:
    double p = 0.5; ///< probability of failure
    double q = 0.5; ///< probability of success
    double logProb = -M_LN2; ///< log(p)
    double log1mProb = -M_LN2; ///< log(q)

private:
    T r = 1; ///< number of failures until the experiment is stopped
    double pdfCoef = -M_LN2; ///< coefficient for faster pdf calculation
    double qDivP = 1; ///< q / p
    static constexpr int tableSize = 16;
    double table[tableSize];
    GammaRand<float> GammaRV{};

protected:
    NegativeBinomialDistribution(T number = 1, double probability = 0.5);

public:
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    IntType MinValue() const override { return 0; }
    IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

protected:
    void SetParameters(T number, double probability);

public:
    inline double GetProbability() const { return p; }
    inline T GetNumber() const { return r; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;

protected:
    enum GENERATOR_ID {
        TABLE,
        EXPONENTIAL,
        GAMMA_POISSON
    };

    /**
     * @fn GetIdOfUsedGenerator
     * If r is small, we use two different generators for two different cases:
     * If p < 0.08 then the tail is too heavy (probability to be in main body is smaller than 0.75),
     * then we return highest integer, smaller than variate from exponential distribution.
     * Otherwise we choose table method
     * @return id of generator
     */
    GENERATOR_ID GetIdOfUsedGenerator() const
    {
        if ((r < 10 || GammaRV.Mean() > 10) && r == std::round(r))
            return (p < 0.08) ? EXPONENTIAL : TABLE;
        return GAMMA_POISSON;
    }

    IntType variateGeometricByTable() const;
    IntType variateGeometricThroughExponential() const;
private:
    IntType variateByTable() const;
    IntType variateThroughExponential() const;
    IntType variateThroughGammaPoisson() const;

public:
    IntType Variate() const override;
    void Sample(std::vector<IntType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
public:
    /**
     * @fn FitProbabilityBayes
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    BetaRand<> FitProbabilityBayes(const std::vector<IntType> &sample, const BetaDistribution<> &priorDistribution, bool MAP = false);
};


/**
 * @brief The NegativeBinomialRand class <BR>
 * Negative binomial distribution
 */
template< typename IntType = int, typename T = double>
class RANDLIBSHARED_EXPORT NegativeBinomialRand : public NegativeBinomialDistribution<IntType, T>
{
public:
    NegativeBinomialRand(T number = 1, double probability = 0.5) : NegativeBinomialDistribution<IntType, T>(number, probability) {}
    String Name() const override;

    using NegativeBinomialDistribution<IntType, T>::SetParameters;

    static constexpr char TOO_SMALL_VARIANCE[] = "Sample variance should be greater than sample mean";

    /**
     * @fn Fit
     * set number and probability, estimated via maximum-likelihood method
     * @param sample
     */
    void Fit(const std::vector<IntType> &sample);
};

template < typename IntType = int >
using PascalRand = NegativeBinomialRand<IntType, int>;

template < typename IntType = int >
using PolyaRand = NegativeBinomialRand<IntType, double>;

#endif // NEGATIVEBINOMIALRAND_H
