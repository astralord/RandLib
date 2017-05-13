#ifndef LIMITINGDISTRIBUTION_H
#define LIMITINGDISTRIBUTION_H

#include "ContinuousDistribution.h"

/**
 * @brief The LimitingDistribution class <BR>
 * Abstract distribution, which has common properties of
 * Stable and Geometric-Stable distributions, due to their similarity
 */
class RANDLIBSHARED_EXPORT LimitingDistribution : public ContinuousDistribution
{
protected:
    double alpha = 2; ///< characteristic exponent α
    double beta = 0; ///< skewness β
    double mu = 0; ///< location μ
    double gamma = M_SQRT2; ///< scale γ
    double alphaInv = M_SQRT1_2; /// 1/α
    double logGamma = -0.5 * M_LN2; ///< log(γ)
    double logGammaPi_2 = M_LNPI - 1.5 * M_LN2; ///< log(γπ/2)

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
