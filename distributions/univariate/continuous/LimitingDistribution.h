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
    double xi, S, alphaInv, zeta; /// coefficients for common α
    double logsigmaPi_2; /// log(σπ/2)

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
    inline double GetInvExponent() const { return alphaInv; }

    double Mean() const override;

protected:
    std::complex<double> psi(double t) const;
};

#endif // LIMITINGDISTRIBUTION_H
