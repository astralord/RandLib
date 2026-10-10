#ifndef NONCENTRALTRAND_H
#define NONCENTRALTRAND_H

#include "ContinuousDistribution.h"
#include "StudentTRand.h"

/**
 * @brief The NoncentralTRand class <BR>
 * Notation: Noncentral-t(ν, μ)
 */
class RANDLIBSHARED_EXPORT NoncentralTRand : public ContinuousDistribution
{
    double nu, mu;
    long double PhiMu, PhimMu;
    double sqrt1p2oNu;
    /// starting point at most weight
    int startingPoint;
    double logHalfMuSq, halfMuSq;
    double lgammaStartingPointpHalf, lgammaStartingPointp1;

    struct nuStruct {
        /// 0.5 ν
        double halfNu;
        double logHalfNu, lgammaHalfNu;
        double lgamma1, lgamma2;
    } nuCoefs, nup2Coefs;

    StudentTRand T;

public:
    explicit NoncentralTRand(double degree = 1, double noncentrality = 0);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double degree, double noncentrality);
    /**
     * @fn GetDegree
     * @return degree ν
     */
    inline double GetDegree() const { return nu; }
    /**
     * @fn GetNoncentrality
     * @return noncentrality μ
     */
    inline double GetNoncentrality() const { return mu; }

private:
    double cdfSeries(const double &x, const nuStruct &degreeCoef, double noncentrality) const;
    double cdfComplSeries(const double &x, const nuStruct &degreeCoef, double noncentrality) const;
    /**
     * @fn logPdfAtZero
     * @return log(f(0))
     */
    double logPdfAtZero() const;
    /**
     * @fn pdfCommon
     * @param x
     * @param noncentrality
     * @return pdf for x ≠ 0
     */
    double pdfCommon(const double & x, double noncentrality) const;
public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};

#endif // NONCENTRALTRAND_H
