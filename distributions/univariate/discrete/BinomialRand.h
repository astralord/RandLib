#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The BinomialDistribution class
 * Abstract class for Binomial distribution
 *
 * Notation: X ~ Bin(n, p)
 *
 * Related distributions:
 * if X ~ Bin(1, p), then X ~ Bernoulli(p)
 */
class RANDLIBSHARED_EXPORT BinomialDistribution : public DiscreteDistribution
{
protected:
    double p, q;
    double logProb, log1mProb;

private:
    int n;
    double np;
    /// log(n!)
    double lgammaNp1;

    double delta1, delta2;
    double sigma1, sigma2, c;
    double a1, a2, a3, a4;
    double coefa3, coefa4;

    /// min(p, q) and [n * min(p, q)] / n respectively
    double minpq, pFloor;
    /// log(pFloor) and log(1 - pFloor)
    double logPFloor, logQFloor;
    /// min(p, q) - pFloor
    double pRes;
    /// [n * p] and [n * q] if p < 0.5, otherwise - vice versa
    double npFloor, nqFloor;
    /// log(P(npFloor))
    double logPnpInv;

    GeometricRand G;

public:
    BinomialDistribution(int number, double probability);
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return n; }

private:
    void SetGeneratorConstants();

protected:
    void SetParameters(int number, double probability);

public:
    inline double GetNumber() const { return n; }
    inline double GetProbability() const { return p; }

private:
    /**
     * @brief logProbFloor
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

    inline double GetLogFactorialN() const { return lgammaNp1; }

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @brief FitProbabilityMLE
     * Fit probability p with maximum-likelihood estimation
     * @param sample
     */
    void FitProbabilityMLE(const std::vector<int> &sample);

    /**
     * @brief FitProbabilityMM
     * Fit probability p, using method of moments
     * @param sample
     */
    void FitProbabilityMM(const std::vector<int> &sample);

    /**
     * @brief FitProbabilityBayes
     * Fit probability p with prior assumption p ~ Beta(α, β)
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    BetaRand FitProbabilityBayes(const std::vector<int> &sample, const BetaDistribution & priorDistribution);

    /**
     * @brief FitProbabilityMinimax
     * Fit probability p with minimax estimator
     * @param sample
     * @return posterior distribution
     */
    BetaRand FitProbabilityMinimax(const std::vector<int> &sample);
};


/**
 * @brief The BinomialRand class
 * Binomial distribution
 */
class RANDLIBSHARED_EXPORT BinomialRand : public BinomialDistribution
{
public:
    BinomialRand(int number, double probability) : BinomialDistribution(number, probability) {}
    std::string Name() const override;
    using BinomialDistribution::SetParameters;
};


#endif // BINOMIALRAND_H
