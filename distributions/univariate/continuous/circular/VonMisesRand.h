#ifndef VONMISESRAND_H
#define VONMISESRAND_H

#include "CircularDistribution.h"

/**
 * @brief The VonMisesRand class <BR>
 * Von-Mises distribution
 *
 * Notation: X ~ Von-Mises(μ, k)
 */
class RANDLIBSHARED_EXPORT VonMisesRand : public CircularDistribution
{
    double k = 1; ///< concentration
    double logI0k = 1; ///< log(I_0(k)), where I_0 stands for modified Bessel function of the first kind of order 0
    double s = M_PI / M_E; ///< generator coefficient
    int p = 12; /// coefficient for faster cdf calculation
    static constexpr double CK = 10.5;

public:
    VonMisesRand(double location, double concentration);

    String Name() const override;

    void SetConcentration(double concentration);
    inline double GetConcentration() const { return k; }

private:
    double cdfSeries(double x) const;
    double cdfErfcAux(double x) const;
    double cdfErfc(double x) const;
    double ccdfErfc(double x) const;

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;

    double CircularMean() const override;
    double CircularVariance() const override;
    double Median() const override;
    double Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // VONMISESRAND_H
