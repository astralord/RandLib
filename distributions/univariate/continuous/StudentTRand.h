#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "GammaRand.h"

/**
 * @brief The StudentTRand class
 * Student's t-distribution
 * X ~ t(v)
 */
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousDistribution
{
    int v;
    double mu, sigma;
    ChiSquaredRand Y;
    double pdfCoef;
    double vp1Half; /// 0.5 * (v + 1)

public:
    explicit StudentTRand(int degree, double location = 0.0, double scale = 1.0);
    std::string name() override;

    void setDegree(int degree);
    void setLocation(double location);
    void setScale(double scale);
    inline int getDegree() const { return v; }
    inline double getLocation() const { return mu; }
    inline double getScale() const { return sigma; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;
    
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // STUDENTTRAND_H
