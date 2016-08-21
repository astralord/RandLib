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

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return -R; }
    double MaxValue() const override { return R; }

    void SetRadius(double radius);
    inline double GetRadius() const { return R; }
    
public:
    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    double Entropy() const;
};

#endif // WIGNERSEMICIRCLERAND_H
