#ifndef SingularDistribution_H
#define SingularDistribution_H

#include "../UnivariateProbabilityDistribution.h"

/**
 * @brief The SingularDistribution class
 */
class RANDLIBSHARED_EXPORT SingularDistribution : public UnivariateProbabilityDistribution<double>
{
public:
    SingularDistribution();
    virtual ~SingularDistribution() {}

    double Hazard(double) const override;
    double Mode() const override;

protected:
    double ExpectedValue(const std::function<double (double)> &, double) const override;

};


#endif // SingularDistribution_H
