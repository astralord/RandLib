#include "ErlangRand.h"

ErlangRand::ErlangRand(int shape, double rate)
{
    setParameters(shape, rate);
}

std::string ErlangRand::name()
{
    return "Erlang(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

void ErlangRand::setParameters(size_t shape, double rate)
{
    GammaRand::setParameters(shape, 1.0 / rate);
}
