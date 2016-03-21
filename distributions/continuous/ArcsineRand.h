#ifndef ARCSINERAND_H
#define ARCSINERAND_H

#include "BetaRand.h"

/**
 * @brief The ArcsineRand class
 */
class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
    double cdfCoef; /// sin(pi * beta) / pi

public:
    ArcsineRand(double minValue = 0, double maxValue = 1, double shape = 0.5);
    std::string name() override;

    void setShape(double shape);
    inline double getShape() const { return beta; }

protected:
    /// prohibit to use beta's getters and setters
    using BetaRand::setShapes;
    using BetaRand::setAlpha;
    using BetaRand::setBeta;
    
public:
    double F(double x) const override;
    double Quantile(double p) const;
    double Mode() const override;
};

#endif // ARCSINERAND_H
