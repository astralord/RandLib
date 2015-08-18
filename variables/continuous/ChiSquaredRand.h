#ifndef CHISQUAREDRAND_H
#define CHISQUAREDRAND_H

#include "ContinuousRand.h"
#include "GammaRand.h"

/**
 * @brief The ChiSquaredRand class
 */
class RANDLIBSHARED_EXPORT ChiSquaredRand : public GammaRand
{
    int k;
    double pdfCoef, cdfCoef;

public:
    ChiSquaredRand(int degree = 1);
    virtual void setName() override;

    void setDegree(int degree);
    inline int getDegree() const { return k; }

    double f(double x) const override;
    double F(double x) const override;

    double E() const override { return k; }
    double Var() const override { return k + k; }

    inline double Mode() const { return std::max(k - 2, 0); }
    inline double Skewness() const { return std::sqrt(8.0 / k); }
    inline double ExcessKurtosis() const { return 12.0 / k; }
};

#endif // CHISQUAREDRAND_H
