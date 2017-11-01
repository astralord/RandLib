#ifndef PARETORAND_H
#define PARETORAND_H

#include "ExponentialRand.h"

/**
 * @brief The ParetoRand class <BR>
 * Pareto distribution
 *
 * Notation: X ~ Pareto(α, σ)
 *
 * Related distributions: <BR>
 * ln(X / σ) ~ Exp(α)
 */
class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousDistribution
{
    double alpha = 1; ///< shape α
    double sigma = 1; ///< scale σ
    double logAlpha = 0; ///< log(α)
    double logSigma = 0; ///< log(σ)

public:
    ParetoRand(double shape = 1, double scale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return sigma; }
    double MaxValue() const override { return INFINITY; }

    void SetShape(double shape);
    void SetScale(double scale);
    inline double GetShape() const { return alpha; }
    inline double GetScale() const { return sigma; }
    inline double GetLogShape() const { return logAlpha; }
    inline double GetLogScale() const { return logSigma; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;

private:
    static double variateForAlphaOne();
    static double variateForAlphaTwo();
    static double variateForGeneralAlpha(double shape);

public:
    static double StandardVariate(double shape);
    static double Variate(double shape, double scale);
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    inline double Entropy() const;

    /**
     * @fn Fit
     * @param sample
     */
    void Fit(const std::vector<double> &sample);
};

#endif // PARETORAND_H
