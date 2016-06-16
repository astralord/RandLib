#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The BinomialRand class
 * Binomial distribution
 * X ~ Bin(n, p)
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

    double pFloor; /// [n * min(p, q)] / n
    double logPFloor, logQFloor; /// log(pFloor) and log(1 - pFloor)
    double pRes; /// min(p, q) - pFloor
    double npFloor, nqFloor; /// [n * p] and [n * q] if p < 0.5, otherwise - vice versa
    double logPnpInv; /// log(P(npFloor))

    static constexpr double generatorEdge = 7.0;

    GeometricRand G;

public:
    BinomialRand(int number, double probability);
    std::string name() const override;

private:
    void setGeneratorConstants();

public:
    void setParameters(int number, double probability);
    inline double getNumber() const { return n; }
    inline double getProbability() const { return p; }

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

private:
    int variateRejection() const;
    int variateWaiting(int number) const;

public:
    int variate() const override;
    static int variate(int n, double p);

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    
    double Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool checkValidity(const std::vector<double> &sample);

    bool fitProbabilityMLE(const std::vector<double> &sample);
    bool fitProbabilityMM(const std::vector<double> &sample);

    /**
     * @brief fitProbability_Bayes
     * fit probability p with assumption p ~ Beta(\alpha, \beta)
     * @param sample
     * @param priorDistribution
     * @return posterior distribution
     */
    bool fitProbabilityBayes(const std::vector<double> &sample, BetaRand &priorDistribution);
};


#endif // BINOMIALRAND_H
