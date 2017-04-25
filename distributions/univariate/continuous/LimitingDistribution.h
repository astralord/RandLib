#ifndef LIMITINGDISTRIBUTION_H
#define LIMITINGDISTRIBUTION_H

#include "ContinuousDistribution.h"

/**
 * @brief The LimitingDistribution class
 * Abstract distribution, which has common properties of
 * Stable and Geometric-stable distributions, due to their similarity
 */
class RANDLIBSHARED_EXPORT LimitingDistribution : public ContinuousDistribution
{
protected:
    double alpha, beta, mu, sigma;
    double alphaInv;

    /// log(σ)
    double logSigma;

    /// log(σπ/2)
    double logsigmaPi_2;

public:
    LimitingDistribution(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~LimitingDistribution() {}

    void SetParameters(double exponent, double skewness);
    void SetLocation(double location);
    void SetScale(double scale);

    inline double GetExponent() const { return alpha; }
    inline double GetSkewness() const { return beta; }
    inline double GetScale() const { return sigma; }
    inline double GetLocation() const { return mu; }
    inline double GetLogScale() const { return logSigma; }

    double Mean() const override;

protected:
    std::complex<double> psi(double t) const;
};

#endif // LIMITINGDISTRIBUTION_H
