#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "LimitingDistribution.h"
#include "StableRand.h"

/**
 * @brief The GeometricStableRand class
 */
class RANDLIBSHARED_EXPORT GeometricStableRand : public LimitingDistribution
{
    StableRand Z;

protected:
    double k; /// asymmetry coefficient for alpha = 2
    double kInv, kSq; /// 1 / k and k * k
    double pdfCoef; /// 1 / (sigma * (k + 1 / k))
    double cdfCoef; /// 1 / (1 + k * k)

public:
    GeometricStableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~GeometricStableRand() {}
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return (alpha < 1 && beta == 1) ? SEMIFINITE_T : INFINITE_T; }
    double MinValue() const override { return (alpha < 1 && beta == 1) ? 0 : -INFINITY; }
    double MaxValue() const override { return INFINITY; }

private:
    void setAsymmetry();
public:
    void setParameters(double exponent, double skewness);
    void setLocation(double location);
    void setScale(double scale);

protected:
    double pdfLaplace(double x) const;
    double cdfLaplace(double x) const;

public:
    double f(double x) const override;
    double F(double x) const override;

private:
    double variateForAlphaEqualOne() const;
    double variateForCommonAlpha() const;
public:
    double variate() const override;
    
    void sample(std::vector<double> &outputData) const override;

    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // GEOMETRICSTABLERAND_H
