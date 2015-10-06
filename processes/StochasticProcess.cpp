#include "StochasticProcess.h"

StochasticProcess::StochasticProcess(double deltaT)
{
    dt = std::max(deltaT, 0.0);
}
