#ifndef IRWINHALLRAND_H
#define IRWINHALLRAND_H

#include "UniformRand.h"

/**
 * @brief The IrwinHallRand class
 */
class RANDLIBSHARED_EXPORT IrwinHallRand : public ContinuousRand
{
    UniformRand U;
    double pdfCoef, cdfCoef;
    int n;
public:
    IrwinHallRand(int number);
    std::string name() override;

    void setNumber(int number);
    inline int getNumber();

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

#endif // IRWINHALLRAND_H
