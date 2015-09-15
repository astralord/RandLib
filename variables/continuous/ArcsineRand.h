#ifndef ARCSINERAND_H
#define ARCSINERAND_H

#include "BetaRand.h"

class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
    double a, b;
    double bma; /// b - a
    double pdfCoef; /// sin(pi * beta) / pi

public:
    ArcsineRand(double minValue, double maxValue, double shape);
    virtual std::string name() override;

    void setSupport(double minValue, double maxValue);
    void setShape(double shape);
    double getMin() { return a; }
    double getMax() { return b; }
    double getShape() { return beta; }

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

    void sample(QVector<double> &outputData);

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const;
    
    double Median() const override;
    double Mode() const override;
};

#endif // ARCSINERAND_H
