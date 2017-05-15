#ifndef PLANCKRAND_H
#define PLANCKRAND_H

#include "GammaRand.h"
#include "../discrete/ZetaRand.h"

/**
 * @brief The PlanckRand class <BR>
 * Planck distribution
 *
 * Notation: X ~ Planck(a, b)
 */
class RANDLIBSHARED_EXPORT PlanckRand : public ContinuousDistribution
{
    double a = 1; ///< shape
    double b = 1; ///< scale
    double pdfCoef = M_LN2 + M_LN3 - 2 * M_LNPI; ///< coefficient for faster pdf calculations

    ZetaRand Z{};
    GammaRand G{2};

public:
    PlanckRand(double shape, double scale);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale);
    inline double GetShape() const { return a; }
    inline double GetScale() const { return b; }

private:
    double leveledPdf(double t) const;

public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // PLANCKRAND_H
