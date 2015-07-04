#include "BasicRandGenerator.h"
#include <time.h>

unsigned long BasicRandGenerator::startPoint = 123456789;

BasicRandGenerator::BasicRandGenerator()
{
    startPoint ^= time(0);
    X = Y = Z = W = startPoint;
    C = 0;
    startPoint = getRand();
}

unsigned long BasicRandGenerator::getRand()
{
    Y ^= Y << 5;
    Y ^= Y >> 7;
    Y ^= Y << 22;

    int t = Z + W + C;
    Z = W;
    C = t < 0;
    W = t & 2147483647;
    X += 1411392427;

    return X + Y + W;
}
