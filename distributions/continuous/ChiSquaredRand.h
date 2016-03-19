#ifndef CHISQUAREDRAND_H
#define CHISQUAREDRAND_H

#include "GammaRand.h"

/**
 * @brief The ChiSquaredRand class
 */
class RANDLIBSHARED_EXPORT ChiSquaredRand : public GammaRand
{
public:
    explicit ChiSquaredRand(int degree = 1);
    std::string name() override;

    void setDegree(int degree);
    inline int getDegree() const;
    
protected:
    /// prohibit to use gamma's public getters and setters
    using GammaRand::setParameters;
    using GammaRand::getShape;
    using GammaRand::getScale;
};

#endif // CHISQUAREDRAND_H
