#ifndef CHIRAND_H
#define CHIRAND_H

#include "ContinuousDistribution.h"
#include "GammaRand.h"

/**
 * @brief The ChiRand class
 */
class RANDLIBSHARED_EXPORT ChiRand : public ContinuousDistribution
{
    int v;
    ChiSquaredRand X;

public:
    explicit ChiRand(int degree);
    void setDegree(int degree);
    inline int getDegree() { return v; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

private:
    double varianceImpl(double mean) const;
    double skewnessImpl(double mean, double sigma) const;

public:
    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // CHIRAND_H
