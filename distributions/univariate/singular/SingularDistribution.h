#ifndef SingularDistribution_H
#define SingularDistribution_H

#include "../UnivariateDistribution.h"

/**
 * @brief The SingularDistribution class <BR>
 * Abstract class for all singular distributions
 */
class RANDLIBSHARED_EXPORT SingularDistribution : public UnivariateDistribution<double>
{
protected:
    SingularDistribution();
    virtual ~SingularDistribution() {}

private:
    double Hazard(double) const override;
    double Mode() const override;
    double ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const override;
    double LikelihoodFunction(const std::vector<double> &sample) const override;
    double LogLikelihoodFunction(const std::vector<double> &sample) const override;
};


#endif // SingularDistribution_H
