#ifndef SKELLAMRAND_H
#define SKELLAMRAND_H

#include "DiscreteRand.h"
#include "PoissonRand.h"

/**
 * @brief The SkellamRand class
 */
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteRand
{
    double mu1, mu2;
    double pmfCoef1; /// exp(-(mu1 + mu2))
    double pmfCoef2; /// sqrt(mu1 / mu2)
    double pmfCoef3; /// 2 * sqrt(mu1 * mu2)

    PoissonRand X, Y;

public:
    SkellamRand(double mean1, double mean2);
    std::string name() override;

    void setMeans(double mean1, double mean2);
    inline double getFirstMean() { return mu1; }
    inline double getSecondMean() { return mu2; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // SKELLAMRAND_H
