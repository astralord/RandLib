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
class RANDLIBSHARED_EXPORT BinomialDistribution : public DiscreteDistribution
{
protected:
    double p = 0.5; ///< probability of success
    double q = 0.5; ///< probability of failure
    double logProb = -M_LN2; ///< log(p)
    double log1mProb = -M_LN2; ///< log(q)

private:
    int n = 1; ///< number of experiments
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

    GeometricRand G{};

protected:
    BinomialDistribution(int number, double probability);

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return n; }

private:
    void SetGeneratorConstants();

protected:
    void SetParameters(int number, double probability);

public:
    inline int GetNumber() const { return n; }
    inline double GetProbability() const { return p; }

private:
    /**
     * @fn logProbFloor
     * @param k
     * @return logarithm of probability to get k if p = pFloor
     */
    double logProbFloor(int k) const;

public:
    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;

private:
    enum GENERATOR_ID {
        BERNOULLI_SUM,
        WAITING,
        REJECTION
    };

    GENERATOR_ID GetIdOfUsedGenerator() const;

    int variateRejection() const;
    int variateWaiting(int number) const;
    static int variateWaiting(int number, double probability);
    static int variateBernoulliSum(int number, double probability);

public:
    int Variate() const override;
    static int Variate(int number, double probability);
    void Sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    int Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

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
    void FitProbability(const std::vector<int> &sample);

    /**
     * @fn FitProbabilityBayes
     * Fit probability p with prior assumption p ~ Beta(α, β)
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    BetaRand FitProbabilityBayes(const std::vector<int> &sample, const BetaDistribution & priorDistribution);

    /**
     * @fn FitProbabilityMinimax
     * Fit probability p with minimax estimator
     * @param sample
     * @return posterior distribution
     */
    BetaRand FitProbabilityMinimax(const std::vector<int> &sample);
};


/**
 * @brief The BinomialRand class <BR>
 * Binomial distribution
 */
class RANDLIBSHARED_EXPORT BinomialRand : public BinomialDistribution
{
public:
    BinomialRand(int number = 1, double probability = 0.5) : BinomialDistribution(number, probability) {}
    String Name() const override;
    using BinomialDistribution::SetParameters;
};


#endif // BINOMIALRAND_H
