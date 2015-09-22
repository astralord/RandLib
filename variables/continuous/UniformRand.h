#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include "ContinuousRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The UniformRand class
 * Continuous uniform distribution: X ~ U(a, b)
 */
class RANDLIBSHARED_EXPORT UniformRand : public ContinuousRand
{
    double a, b;
    double c; /// 1 / (b - a)

public:
    UniformRand(double minValue = 0, double maxValue = 1);
    std::string name() override;

    void setBoundaries(double minValue, double maxValue);    
    inline double getMinValue() const { return a; }
    inline double getMaxValue() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double minValue, double maxValue);
    static double standardVariate();

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() const;

    bool fitToData(const QVector<double> &sample);
};

#endif // UNIFORMRAND_H
