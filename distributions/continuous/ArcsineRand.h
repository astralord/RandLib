#ifndef ARCSINERAND_H
#define ARCSINERAND_H

#include "BetaRand.h"

class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
    double a, b;
    double bma; /// b - a
    double pdfCoef; /// sin(pi * beta) / pi

public:
    ArcsineRand(double minValue = 0, double maxValue = 1, double shape = 0.5);
    std::string name() override;

    void setSupport(double minValue, double maxValue);
    void setShape(double shape);
    inline double getMin() const { return a; }
    inline double getMax() const { return b; }
    inline double getShape() const { return beta; }

protected:
    /// prohibit to use beta's getters and setters
    using BetaRand::setParameters;
    using BetaRand::setAlpha;
    using BetaRand::setBeta;
    using BetaRand::getAlpha;
    using BetaRand::getBeta;
    
public:
    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const;
    
    double Median() const override;
    double Mode() const override;
};

#endif // ARCSINERAND_H
