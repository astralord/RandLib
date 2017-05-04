#ifndef VONMISESRAND_H
#define VONMISESRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The VonMisesRand class <BR>
 * Von-Mises distribution
 *
 * Notation: X ~ Von-Mises(μ, k)
 */
class RANDLIBSHARED_EXPORT VonMisesRand : public ContinuousDistribution
{
    double mu, k;
    /// log(I_0(k))
    double logI0k;
     /// generator coefficient
    double s;
    /// cdf coefficient
    int p;
    static constexpr double CK = 10.5;

public:
    VonMisesRand(double location, double concentration);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return mu - M_PI; }
    double MaxValue() const override { return mu + M_PI; }

    void SetLocation(double location);
    void SetConcentration(double concentration);
    inline double GetLocation() const { return mu; }
    inline double GetConcentration() const { return k; }

private:
    double cdfSeries(double x) const;
    double cdfErfc(double x) const;

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // VONMISESRAND_H
