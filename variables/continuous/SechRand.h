#ifndef SECHRAND_H
#define SECHRAND_H

#include "ContinuousRand.h"

class RANDLIBSHARED_EXPORT SechRand : public ContinuousRand
{
public:
    SechRand();
    std::string name() override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double Kurtosis() const override;

    double Entropy() const;
};

#endif // SECHRAND_H
