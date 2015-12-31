#ifndef BINOMIALRAND_H
#define BINOMIALRAND_H

#include "DiscreteRand.h"
#include "GeometricRand.h"

/**
 * @brief The BinomialRand class
 */
class RANDLIBSHARED_EXPORT BinomialRand : public DiscreteRand
{
    int n;
    double p, q;
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
    std::string name() override;

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
    double F(double x) const override;

private:
    double variateRejection() const;
    double variateWaiting(int number) const;

public:
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


#endif // BINOMIALRAND_H
