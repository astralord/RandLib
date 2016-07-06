#ifndef LIMITINGDISTRIBUTION_H
#define LIMITINGDISTRIBUTION_H

#include "ContinuousDistribution.h"

/**
 * @brief The LimitingDistribution class
 */
class RANDLIBSHARED_EXPORT LimitingDistribution : public ContinuousDistribution
{
protected:
    double alpha, beta, mu, sigma;
    double B, S, alphaInv, zeta; /// coefficients for alpha != 1
    double logSigma; /// coefficients for alpha == 1

public:
    LimitingDistribution(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~LimitingDistribution() {}

    void setParameters(double exponent, double skewness);
    void setLocation(double location);
    void setScale(double scale);

    inline double getExponent() const { return alpha; }
    inline double getSkewness() const { return beta; }
    inline double getScale() const { return sigma; }
    inline double getLocation() const { return mu; }

    double Mean() const override;

protected:
    std::complex<double> psi(double t) const;
};

#endif // LIMITINGDISTRIBUTION_H
