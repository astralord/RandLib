#ifndef RAABGREENRAND_H
#define RAABGREENRAND_H

#include "ContinuousRand.h"

class RANDLIBSHARED_EXPORT RaabGreenRand : public ContinuousRand
{
public:
    RaabGreenRand();
    std::string name() override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;

    // TODO: implement all!
};

#endif // RAABGREENRAND_H
