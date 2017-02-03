#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The BinomialRand class
 * Binomial distribution
 *
 * Notation: X ~ Bin(n, p)
 *
 * Related distributions:
 * if X ~ Bin(1, p), then X ~ Bernoulli(p)
 */
class RANDLIBSHARED_EXPORT BinomialRand : public DiscreteDistribution
{
protected:
    double p, q;

private:
    int n;
    double np;

    double delta1, delta2;
    double sigma1, sigma2, c;
    double a1, a2, a3, a4;
    double coefa3, coefa4;

    double minpq, pFloor; /// min(p, q) and [n * min(p, q)] / n respectively
    double logPFloor, logQFloor; /// log(pFloor) and log(1 - pFloor)
    double pRes; /// min(p, q) - pFloor
    double npFloor, nqFloor; /// [n * p] and [n * q] if p < 0.5, otherwise - vice versa
    double logPnpInv; /// log(P(npFloor))

    GeometricRand G;

public:
    BinomialRand(int number, double probability);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return n; }

private:
    void SetGeneratorConstants();

public:
    void SetParameters(int number, double probability);
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
    double P(int k) const override;
    double F(int k) const override;
    double S(int k) const override;

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
    double Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @brief FitProbabilityMLE
     * Fit probability p with maximum-likelihood estimation
     * @param sample
     * @return
     */
    bool FitProbabilityMLE(const std::vector<int> &sample);

    /**
     * @brief FitProbabilityMM
     * Fit probability p, using method of moments
     * @param sample
     * @return
     */
    bool FitProbabilityMM(const std::vector<int> &sample);

    /**
     * @brief FitProbability_Bayes
     * Fit probability p with prior assumption p ~ Beta(α, β)
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    bool FitProbabilityBayes(const std::vector<int> &sample, BetaRand &priorDistribution);

    /**
     * @brief FitProbabilityMinimax
     * Fit probability p with minimax estimator
     * @param sample
     * @return
     */
    bool FitProbabilityMinimax(const std::vector<int> &sample);
};


#endif // BINOMIALRAND_H
