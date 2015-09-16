#ifndef PARETORAND_H
#define PARETORAND_H

#include "ExponentialRand.h"

/**
 * @brief The ParetoRand class
 */
class RANDLIBSHARED_EXPORT ParetoRand : public ContinuousRand
{
    double xm, alpha;
    double alphaInv;
    double pdfCoef;

public:
    ParetoRand(double shape, double scale);
    std::string name() override;

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    inline double getShape() const { return xm; }
    inline double getScale() const { return alpha; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double shape, double scale);

private:
    static double variateForAlphaOne();
    static double variateForAlphaTwo();
    static double variateForCommonAlpha(double shape);

public:
    void sample(QVector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() const;

    bool fitToData(const QVector<double> &sample);
};

#endif // PARETORAND_H
