#ifndef PARETORAND_H
#define PARETORAND_H

#include "ExponentialRand.h"

/**
 * @brief The ParetoRand class
 */
class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousDistribution
{
    double xm, alpha;
    double alphaInv;
    double pdfCoef;

public:
    ParetoRand(double shape = 1, double scale = 1);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return xm; }
    double MaxValue() const override { return INFINITY; }

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    inline double getShape() const { return alpha; }
    inline double getScale() const { return xm; }

    double f(double x) const override;
    double F(double x) const override;

private:
    static double variateForAlphaOne();
    static double variateForAlphaTwo();
    static double variateForCommonAlpha(double shape);

public:
    static double standardVariate(double shape);
    static double variate(double shape, double scale);
    double variate() const override;

    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() const;

    bool fitShapeAndScaleMLE(const std::vector<double> &sample);
};

#endif // PARETORAND_H
