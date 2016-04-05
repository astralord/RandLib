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

protected:
    double sigma;
    double sigmaSqInv; // 1.0 / sigma^2

public:
    explicit ChiRand(int degree, double scale = 1.0);
    std::string name() override;

    void setDegree(int degree);
    void setScale(double scale);
    inline int getDegree() { return X.getDegree(); }
    inline double getScale() { return sigma; }

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


/**
 * @brief The MaxwellBoltzmannRand class
 * Maxwell-Boltzmann distribution
 * X ~ MB(a)
 *
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public ChiRand
{
public:
    explicit MaxwellBoltzmannRand(double scale);
    std::string name() override;

    double f(double x) const override;
    double F(double x) const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The RayleighRand class
 */
class RANDLIBSHARED_EXPORT RayleighRand : public ChiRand
{
public:
    explicit RayleighRand(double scale = 1);
    std::string name() override;

    double f(double x) const override;
    double F(double x) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool checkValidity(const std::vector<double> &sample);

    bool fitScaleMLE(const std::vector<double> &sample);
    bool fitScaleUMVU(const std::vector<double> &sample);
};

#endif // CHIRAND_H
