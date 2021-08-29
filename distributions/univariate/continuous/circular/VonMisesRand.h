#ifndef VONMISESRAND_H
#define VONMISESRAND_H

#include "CircularDistribution.h"

/**
 * @brief The VonMisesRand class <BR>
 * Von-Mises distribution
 *
 * Notation: X ~ Von-Mises(μ, k)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT VonMisesRand : public CircularDistribution<RealType>
{
    double k = 1; ///< concentration
    double logI0k = 1; ///< log(I_0(k)), where I_0 stands for modified Bessel function of the first kind of order 0
    double s = M_PI / M_E; ///< generator coefficient
    int p = 12; /// coefficient for faster cdf calculation
    static constexpr double CK = 10.5;

public:
    VonMisesRand(double location = 0, double concentration = 1);

    String Name() const override;

    void SetConcentration(double concentration);
    inline double GetConcentration() const { return k; }

private:
    double cdfSeries(double x) const;
    double cdfErfcAux(RealType x) const;
    double cdfErfc(RealType x) const;
    double ccdfErfc(RealType x) const;

public:
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    long double CircularMean() const override;
    long double CircularVariance() const override;
    RealType Median() const override;
    RealType Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // VONMISESRAND_H
