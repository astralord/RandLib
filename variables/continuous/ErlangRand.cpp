#include "ErlangRand.h"

ErlangRand::ErlangRand(int shape, double rate)
{
    setParameters(shape, rate);
    name();
}

std::string ErlangRand::name()
{
    return "Erlang(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

void ErlangRand::setParameters(int shape, double rate)
{
    GammaRand::setParameters(static_cast<double>(shape), 1.0 / rate);
}
