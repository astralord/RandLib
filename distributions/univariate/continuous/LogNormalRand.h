#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"

/**
 * @brief The LogNormalRand class
 * Log-normal distribution
 */
class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousDistribution
{
    NormalRand X;

public:
    LogNormalRand(double location = 0, double squaredScale = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return X.Mean(); }
    inline double GetScale() const { return X.GetScale(); }

    double f(double x) const override;
    double F(double x) const override;

    static double StandardVariate();
    static double Variate(double location, double scale);
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    double logAverage(const std::vector<double> &sample);
    double logSecondMoment(const std::vector<double> &sample);

public:
    /// Method of moments
    bool FitLocationMM(const std::vector<double> &sample);
    bool FitScaleMM(const std::vector<double> &sample);
    bool FitMM(const std::vector<double> &sample);

    /// Maximum-likelihod estimation
    bool FitLocationMLE(const std::vector<double> &sample);
    bool FitScaleMLE(const std::vector<double> &sample);
    bool FitMLE(const std::vector<double> &sample);

    /// Bayesian estimation
    bool FitLocationBayes(const std::vector<double> &sample, NormalRand &priorDistribution);
    bool FitScaleBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution);
    bool FitBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution);
};

#endif // LOGNORMALRAND_H
