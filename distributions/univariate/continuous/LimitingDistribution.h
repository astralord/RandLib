#ifndef LIMITINGDISTRIBUTION_H
#define LIMITINGDISTRIBUTION_H

#include "ContinuousDistribution.h"

/**
 * @brief The LimitingDistribution class
 * Abstract distribution, which has common properties of
 * Stable and Geometric-Stable distributions, due to their similarity
 */
class RANDLIBSHARED_EXPORT LimitingDistribution : public ContinuousDistribution
{
protected:
    double alpha, beta, mu, gamma;
    double alphaInv;

    /// log(γ)
    double logGamma;

    /// log(γπ/2)
    double logGammaPi_2;

public:
    LimitingDistribution(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~LimitingDistribution() {}
protected:
    void SetParameters(double exponent, double skewness);
public:
    void SetLocation(double location);
    void SetScale(double scale);

    inline double GetExponent() const { return alpha; }
    inline double GetSkewness() const { return beta; }
    inline double GetScale() const { return gamma; }
    inline double GetLocation() const { return mu; }
    inline double GetLogScale() const { return logGamma; }

    double Mean() const override;

protected:
    std::complex<double> psi(double t) const;
};

#endif // LIMITINGDISTRIBUTION_H
