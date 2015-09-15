#ifndef ERLANGRAND_H
#define ERLANGRAND_H

#include "GammaRand.h"

/**
 * @brief The ErlangRand class
 * Erlang distibution:
 * X ~ Erlang(k, l)
 * X ~ Y_1 + Y_2 + ... + Y_k, where Y_i ~ Exp(l)
 * X ~ Gamma(k, 1/l)
 */
class RANDLIBSHARED_EXPORT ErlangRand : public GammaRand
{
public:
    ErlangRand(int shape = 1, double rate = 1);
    virtual std::string name() override;

    void setParameters(int shape, double rate);
    inline int getShape() const { return static_cast<int>(GammaRand::getShape()); }
    inline double getRate() const { return 1.0 / GammaRand::getScale(); }
    
protected:
    /// prohibit to use gamma's public getters and setters
    using GammaRand::setParameters;
    using GammaRand::getShape;
    using GammaRand::getScale;
};

#endif // ERLANGRAND_H
