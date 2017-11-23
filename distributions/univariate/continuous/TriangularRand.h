#ifndef TRIANGULARRAND_H
#define TRIANGULARRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The TriangularRand class <BR>
 * Triangular distribution
 *
 * Notation: X ~ Tri(a, b, c)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT TriangularRand : public ContinuousDistribution<RealType>
{
    double a = 0; ///< min value
    double b = 2; ///< max value
    double c = 1; ///< mode
    double constForGenerator = 1; ///< (c - a) / (b - a)
    double coefGenerator1 = 1; ///< (b - a) * (c - a)
    double coefGenerator2 = 1; ///< (b - a) * (b - c)

    void SetConstantsForGenerator();

public:
    TriangularRand(double lowerLimit = 0, double mode = 0.5, double upperLimit = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return a; }
    RealType MaxValue() const override { return b; }

    void SetParameters(double lowerLimit, double mode, double upperLimit);

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

public:
    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // TRIANGULARRAND_H
