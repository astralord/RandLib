#ifndef SINGULARRAND_H
#define SINGULARRAND_H

#include "../RandomVariable.h"

/**
 * @brief The SingularRand class
 */
class RANDLIBSHARED_EXPORT SingularRand : public RandomVariable
{
public:
    SingularRand();
    virtual ~SingularRand() {}

    double Quantile(double p) const override;

protected:
    double Hazard(double) const override;
    double ExpectedValue(const std::function<double (double)> &, double) const override;

    double Mode() const override;
};


#endif // SINGULARRAND_H
