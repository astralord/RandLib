#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The BinomialDistribution class <BR>
 * Abstract class for Binomial distribution
 *
 * Notation: X ~ Bin(n, p)
 *
 * Related distributions: <BR>
 * If X ~ Bin(1, p), then X ~ Bernoulli(p) <BR>
 * X ~ Multin(n, 1 - p, p)
 */
template< typename IntType = int >
class RANDLIBSHARED_EXPORT BinomialDistribution : public DiscreteExponentialFamily<IntType, double>
{
protected:
    double p = 0.5; ///< probability of success
    double q = 0.5; ///< probability of failure
    double logProb = -M_LN2; ///< log(p)
    double log1mProb = -M_LN2; ///< log(q)

private:
    IntType n = 1; ///< number of experiments
    double np = 0.5; ///< n * p
    double lfactn = 0; ///< log(n!)

    double delta1{}, delta2{};
    double sigma1{}, sigma2{}, c{};
    double a1{}, a2{}, a3{}, a4{};
    double coefa3{}, coefa4{};

    double minpq = 0.5; ///< min(p, q)
    double pFloor = 0; ///< [n * min(p, q)] / n
    double logPFloor = -INFINITY; ///< log(pFloor)
    double logQFloor = 0; ///< log(1 - pFloor)
    double pRes = 0.5; ///< min(p, q) - pFloor
    double npFloor = 0; ///< [n * min(p, q)]
    double nqFloor = 0; ///< [n * max(p, q)]
    double logPnpInv = 0; ///< log(P([npFloor)) if p = pFloor

    GeometricRand<IntType> G{};

protected:
    BinomialDistribution(IntType number, double probability);

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return 0; }
    IntType MaxValue() const override { return n; }

private:
    void SetGeneratorConstants();

protected:
    void SetParameters(IntType number, double probability);

public:
    inline IntType GetNumber() const { return n; }
    inline double GetProbability() const { return p; }

    double SufficientStatistic(IntType x) const override;
    double SourceParameters() const override;
    double SourceToNatural(double sourceParameters) const override;
    double NaturalParameters() const override;
    double LogNormalizer(double theta) const override;
    double LogNormalizerGradient(double theta) const override;
    double CarrierMeasure(IntType x) const override;

private:
    /**
     * @fn logProbFloor
     * @param k
     * @return logarithm of probability to get k if p = pFloor
     */
    double logProbFloor(int k) const;

public:
    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;

private:
    enum GENERATOR_ID {
        BERNOULLI_SUM,
        WAITING,
        REJECTION,
        POISSON
    };

    GENERATOR_ID GetIdOfUsedGenerator() const
    {
        /// if (n is tiny and minpq is big) or p = 0.5 and n is not so large,
        /// we just sum Bernoulli random variables
        if ((n <= 3) || (n <= 13 && minpq > 0.025 * (n + 6)) || (n <= 200 && RandMath::areClose(p, 0.5)))
            return BERNOULLI_SUM;

        /// for small [np] we use simple waiting algorithm
        if ((npFloor <= 12) ||
            (pRes > 0 && npFloor <= 16))
            return WAITING;

        /// otherwise
        return REJECTION;
    }

    IntType variateRejection() const;
    IntType variateWaiting(IntType number) const;
    static IntType variateWaiting(IntType number, double probability, RandGenerator &randGenerator);
    static IntType variateBernoulliSum(IntType number, double probability, RandGenerator &randGenerator);

public:
    IntType Variate() const override;
    static IntType Variate(IntType number, double probability, RandGenerator &randGenerator = ProbabilityDistribution<IntType>::staticRandGenerator);
    void Sample(std::vector<IntType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Median() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

    /**
     * @fn GetLogFactorialN
     * @return log(n!)
     */
    inline double GetLogFactorialN() const { return lfactn; }
    /**
     * @fn GetLogProbability
     * @return log(p)
     */
    inline double GetLogProbability() const { return logProb; }
    /**
     * @fn GetLog1mProbability
     * @return log(1-p)
     */
    inline double GetLog1mProbability() const { return log1mProb; }

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @fn FitProbability
     * Fit probability p with maximum-likelihood estimation
     * @param sample
     */
    void FitProbability(const std::vector<IntType> &sample);

    /**
     * @fn FitProbabilityBayes
     * Fit probability p with prior assumption p ~ Beta(α, β)
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    BetaRand<> FitProbabilityBayes(const std::vector<IntType> &sample, const BetaDistribution<> & priorDistribution, bool MAP = false);

    /**
     * @fn FitProbabilityMinimax
     * Fit probability p with minimax estimator
     * @param sample
     * @return posterior distribution
     */
    BetaRand<> FitProbabilityMinimax(const std::vector<IntType> &sample);
};


/**
 * @brief The BinomialRand class <BR>
 * Binomial distribution
 */
template< typename IntType = int >
class RANDLIBSHARED_EXPORT BinomialRand : public BinomialDistribution<IntType>
{
public:
    BinomialRand(int number = 1, double probability = 0.5) : BinomialDistribution<IntType>(number, probability) {}
    String Name() const override;
    using BinomialDistribution<IntType>::SetParameters;
};


#endif // BINOMIALRAND_H
