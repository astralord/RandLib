#ifndef SingularDistribution_H
#define SingularDistribution_H

#include "../ProbabilityDistribution.h"

/**
 * @brief The SingularDistribution class
 */
class RANDLIBSHARED_EXPORT SingularDistribution : public ProbabilityDistribution
{
public:
    SingularDistribution();
    virtual ~SingularDistribution() {}

protected:
    double Hazard(double) const override;
    double ExpectedValue(const std::function<double (double)> &, double) const override;

    double Mode() const override;
};


#endif // SingularDistribution_H
