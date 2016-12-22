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
    /// parameters for α = 2
    double k; /// asymmetry coefficient
    double kInv, kSq; /// 1 / k and k * k
    double pdfCoef; /// 1 / (σ * (k + 1 / k))
    double cdfCoef; /// 1 / (1 + k * k)

public:
    GeometricStableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~GeometricStableRand() {}

    std::string Name() const override;

    SUPPORT_TYPE SupportType() const override {
        if (alpha < 1) {
            if (beta == 1 && mu >= 0)
                return RIGHTSEMIFINITE_T;
            if (beta == -1 && mu <= 0)
                return LEFTSEMIFINITE_T;
        }
        return INFINITE_T;
    }

    double MinValue() const override {
        if (alpha < 1 && beta == 1 && mu >= 0)
            return 0;
        return -INFINITY;
    }

    double MaxValue() const override {
        if (alpha < 1 && beta == -1 && mu <= 0)
            return 0;
        return INFINITY;
    }

public:
    void SetParameters(double exponent, double skewness, double scale, double location);

protected:
    double pdfLaplace(double x) const;
    double cdfLaplace(double x) const;

private:
    double pdfByLevy(double x) const;
    double pdfByCauchy(double x) const;

public:
    double f(double x) const override;
    double F(double x) const override;

private:
    double variateForUnityExponent() const;
    double variateForCommonExponent() const;
    double variateByLevy(bool positive) const;
    double variateByCauchy() const;
public:
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;
};

#endif // GEOMETRICSTABLERAND_H
