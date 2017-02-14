#ifndef RAISEDCOSINERAND_H
#define RAISEDCOSINERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The RaisedCosineRand class
 *
 * X ~ Raised-cosine(μ, s)
 */
class RANDLIBSHARED_EXPORT RaisedCosineRand : public ContinuousDistribution
{
    double mu, s;
    double s_pi;
public:
    RaisedCosineRand(double location, double scale);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return mu - s; }
    double MaxValue() const override { return mu + s; }

    void SetLocation(double location);
    inline double GetLocation() const { return mu; }
    void SetScale(double scale);
    inline double GetScale() const { return s; }

    double f(double x) const override;
    double F(double x) const override;

    static double StandardVariate();
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};


/**
 * @brief The RaabGreenRand class
 *
 * X ~ Raab-Green()
 *
 * Related distributions:
 * X ~ Raised-cosine(0.0, π)
 */
class RANDLIBSHARED_EXPORT RaabGreenRand : public RaisedCosineRand
{
public:
    RaabGreenRand() : RaisedCosineRand(0.0, M_PI) {}
    std::string Name() const override;
};


#endif // RAISEDCOSINERAND_H
