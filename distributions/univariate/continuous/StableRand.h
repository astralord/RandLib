#ifndef STABLERAND_H
#define STABLERAND_H

#include <functional>
#include "LimitingDistribution.h"

/**
 * @brief The StableRand class
 */
class RANDLIBSHARED_EXPORT StableRand : public LimitingDistribution
{
    double alpham1Inv, alpha_alpham1; /// 1 / (alpha - 1) and alpha / (alpha - 1)
    double xi, integrandCoef;

protected:
    double pdfCoef;

public:
    StableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~StableRand() {}

    std::string name() const override;
    SUPPORT_TYPE supportType() const override {
        if (alpha < 1) {
            if (beta == 1)
                return RIGHTSEMIFINITE_T;
            if (beta == -1)
                return LEFTSEMIFINITE_T;
        }
        return INFINITE_T;
    }
    double MinValue() const override { return (alpha < 1 && beta == 1) ? mu : -INFINITY; }
    double MaxValue() const override { return (alpha < 1 && beta == -1) ? mu : INFINITY; }

    void setParameters(double exponent, double skewness);
    void setScale(double scale);
    
protected:
    /// Probability distribution functions
    double pdfNormal(double x) const;
    double pdfCauchy(double x) const;
    double pdfLevy(double x) const;
private:
    double integrandAuxForAlphaEqualOne(double theta, double xAdj) const;
    double integrandForAlphaEqualOne(double theta, double xAdj) const;
    double pdfForAlphaEqualOne(double x) const;

    double integrandAuxForCommonAlpha(double theta, double xAdj, double xiAdj) const;
    double integrandForCommonAlpha(double theta, double xAdj, double xiAdj) const;
    double pdfForCommonAlpha(double x) const;
public:    
    double f(double x) const override;
    
protected:
    /// Cumulative distribution functions
    double cdfNormal(double x) const;
    double cdfCauchy(double x) const;
    double cdfLevy(double x) const;
private:
    double cdfForCommonAlpha(double x) const;
    double cdfForAlphaEqualOne(double x) const;
public:
    double F(double x) const override;
    
private:
    /// variate
    double variateForCommonAlpha() const;
    double variateForAlphaEqualOne() const;
public:
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

public:
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The HoltsmarkRand class
 * X ~ Stable(1.5, 0.0, scale, location)
 */
class RANDLIBSHARED_EXPORT HoltsmarkRand : public StableRand
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableRand(1.5, 0.0, scale, location) {}
    std::string name() const override;
};

/**
 * @brief The LandauRand class
 * X ~ Stable(1.0, 1.0, scale, location)
 */
class RANDLIBSHARED_EXPORT LandauRand : public StableRand
{
public:
    LandauRand(double scale = 1, double location = 0) : StableRand(1.0, 1.0, scale, location) {}
    std::string name() const override;
};

#endif // STABLERAND_H
