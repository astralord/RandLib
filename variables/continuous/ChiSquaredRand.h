#ifndef CHISQUAREDRAND_H
#define CHISQUAREDRAND_H

#include "ContinuousRand.h"
#include "GammaRand.h"

/**
 * @brief The ChiSquaredRand class
 */
class RANDLIBSHARED_EXPORT ChiSquaredRand : public GammaRand
{
public:
    explicit ChiSquaredRand(int degree = 1);
    virtual std::string name() override;

    void setDegree(int degree);
    inline int getDegree() const { return k + k; }
};

#endif // CHISQUAREDRAND_H
