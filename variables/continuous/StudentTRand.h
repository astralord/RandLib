#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "ContinuousRand.h"
#include "ChiSquaredRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"

/**
 * @brief The StudentTRand class
 */
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousRand
{
    int v;
    ChiSquaredRand Y;
    double pdfCoef;
public:
    explicit StudentTRand(int degree);
    virtual std::string name() override;

    void setDegree(int degree);
    inline int getDegree() const { return v; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(QVector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // STUDENTTRAND_H
