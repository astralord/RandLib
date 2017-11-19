#ifndef FISHERFRAND_H
#define FISHERFRAND_H

#include "BetaPrimeRand.h"

/**
 * @brief The FisherFRand class <BR>
 * F-distribution
 *
 * Notation: X ~ F(d1, d2)
 *
 * Related distributions: <BR>
 * d1/d2 * X ~ B'(d1/2, d2/2)
 */
template < typename RealType = long double >
class RANDLIBSHARED_EXPORT FisherFRand : public ContinuousDistribution<RealType>
{
    int d1 = 2; ///< first degree
    int d2 = 2; ///< second degree
    double a = 0; ///< d1 / 2 - 1;
    double d1_d2 = 1; ///< d1 / d2
    double c = -2; ///< -(d1 + d2) / 2;
    double d2_d1 = 1; ///< d2 / d1
    double pdfCoef = 0; /// < (a + 1) * log(d1/d2) - log(B(d1/2, d2/2))

    BetaPrimeRand<RealType> B{};

public:
    FisherFRand(int degree1, int degree2);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

    void SetDegrees(int degree1, int degree2);
    inline int GetFirstDegree() const { return d1; }
    inline int GetSecondDegree() const { return d2; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // FISHERFRAND_H
