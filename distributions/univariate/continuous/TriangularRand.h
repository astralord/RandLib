#ifndef TRIANGULARRAND_H
#define TRIANGULARRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The TriangularRand class
 * Triangular distribution
 *
 * Notation: X ~ Tri(a, b, c)
 */
class RANDLIBSHARED_EXPORT TriangularRand : public ContinuousDistribution
{
    double a, b, c;
    double constForGenerator; /// (c - a) / (b - a)
    double coefGenerator1; /// (b - a) * (c - a)
    double coefGenerator2; /// (b - a) * (b - c)
    void setConstantsForGenerator();

public:
    TriangularRand(double lowerLimit = 0, double mode = 0.5, double upperLimit = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return b; }

    void setParameters(double lowerLimit, double mode, double upperLimit);
    inline double getLowerLimit() const { return a; }
    inline double getMode() const { return c; }
    inline double getUpperLimit() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

public:
    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // TRIANGULARRAND_H
