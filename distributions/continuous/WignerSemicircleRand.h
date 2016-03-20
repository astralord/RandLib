#ifndef WIGNERSEMICIRCLERAND_H
#define WIGNERSEMICIRCLERAND_H

#include "BetaRand.h"

/**
 * @brief The WignerSemicircleRand class
 */
class RANDLIBSHARED_EXPORT WignerSemicircleRand : public ContinuousDistribution
{
    double R, RSq;
    BetaRand X; /// for generator
    
public:
    explicit WignerSemicircleRand(double radius);
    std::string name() override;

    void setRadius(double radius);
    inline double getRadius() const { return R; }
    
public:
    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    double Entropy() const;
};

#endif // WIGNERSEMICIRCLERAND_H
