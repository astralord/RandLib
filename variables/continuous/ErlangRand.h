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
    ErlangRand(int shape, double rate);

    void setParameters(int shape, double rate);
    inline int getShape() { return static_cast<int>(GammaRand::getShape()); }
    inline double getRate() { return 1.0 / GammaRand::getScale(); }
};

#endif // ERLANGRAND_H
