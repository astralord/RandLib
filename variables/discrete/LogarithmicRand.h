#ifndef LOGARITHMICRAND_H
#define LOGARITHMICRAND_H

#include "DiscreteRand.h"

/**
 * @brief The LogarithmicRand class
 */
class RANDLIBSHARED_EXPORT LogarithmicRand : public DiscreteRand
{
    double p, q;
    double logQInv; /// 1 / log(q)
public:
    LogarithmicRand(double probability);
    std::string name() override;

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Mode() const override;
};

#endif // LOGARITHMICRAND_H
