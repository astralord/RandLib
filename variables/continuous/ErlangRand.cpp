#include "ErlangRand.h"

ErlangRand::ErlangRand(int shape, double rate)
{
    setParameters(shape, rate);
}

std::string ErlangRand::name()
{
    return "Erlang(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

void ErlangRand::setParameters(int shape, double rate)
{
    GammaRand::setParameters(std::max(shape, 1), 1.0 / rate);
}
