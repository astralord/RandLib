#ifndef NAKAGAMIRAND_H
#define NAKAGAMIRAND_H

#include "GammaRand.h"

/**
 * @brief The NakagamiRand class
 */
class RANDLIBSHARED_EXPORT NakagamiRand : public ContinuousDistribution
{
    double m, w;
    double sigma; /// w / m

    GammaRand Y;

public:
    NakagamiRand(double shape = 0.5, double spread = 1);
    std::string name() const override;

    void setParameters(double shape, double spread);
    inline double getShape() const { return m; }
    inline double getSpread() const { return w; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
};


/**
 * @brief The ChiRand class
 */
class RANDLIBSHARED_EXPORT ChiRand : public NakagamiRand
{
protected:
    double sigma;
    double sigmaSqInv; /// 1.0 / sigma^2

public:
    explicit ChiRand(int degree, double scale = 1.0);
    std::string name() const override;

private:
    using NakagamiRand::setParameters;
    using NakagamiRand::getShape;
    using NakagamiRand::getSpread;

public:
    void setParameters(int degree, double scale = 1.0);
    inline int getDegree() const { return 2 * NakagamiRand::getShape(); }
    inline double getScale() const { return sigma; }

private:
    double skewnessImpl(double mean, double sigma) const;

public:
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
    std::string name() const override;

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
    std::string name() const override;

    void setScale(double scale);

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


#endif // NAKAGAMIRAND_H
