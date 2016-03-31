#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"

/**
 * @brief The LogNormalRand class
 */
class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousDistribution
{
    double expMu, expVar;
    NormalRand X;

public:
    LogNormalRand(double location = 0, double scale = 1);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return X.Mean(); }
    inline double getScale() const { return X.getScale(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double logAverage(const std::vector<double> &sample);
    double logSecondMoment(const std::vector<double> &sample);

public:
    bool checkValidity(const std::vector<double> &sample);

    /// Method of moments
    bool fitLocationMM(const std::vector<double> &sample);
    bool fitScaleMM(const std::vector<double> &sample);
    bool fitLocationAndScaleMM(const std::vector<double> &sample);

    /// Maximum-likelihod estimation
    bool fitLocationMLE(const std::vector<double> &sample);
    bool fitScaleMLE(const std::vector<double> &sample);
    bool fitLocationAndScaleMLE(const std::vector<double> &sample);

    /// Bayesian estimation
    bool fitLocationBayes(const std::vector<double> &sample, NormalRand &priorDistribution);
    bool fitScaleBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution);
    bool fitLocationAndScaleBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution);
};

#endif // LOGNORMALRAND_H
