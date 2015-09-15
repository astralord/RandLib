#ifndef ARCSINERAND_H
#define ARCSINERAND_H

#include "BetaRand.h"

class RANDLIBSHARED_EXPORT ArcsineRand : public BetaRand
{
    double a, b;
    double pdfCoef; /// sin(pi * shape) / pi

public:
    ArcsineRand(double minValue, double maxValue, double shape);
    virtual std::string name() override;

    void setSupport(double minValue, double maxValue);
    void setShape(double shape);
    double getMin() { return a; }
    double getMax() { return b; }
    double getShape() { return BetaRand::getBeta(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // ARCSINERAND_H
