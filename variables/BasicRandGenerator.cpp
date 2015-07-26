#include "BasicRandGenerator.h"
#include <time.h>

unsigned long BasicRandGenerator::startPoint = 123456789 ^ time(0);
unsigned long BasicRandGenerator::X = startPoint;
unsigned long BasicRandGenerator::Y = startPoint;
unsigned long BasicRandGenerator::Z = startPoint;
unsigned long BasicRandGenerator::W = startPoint;
bool BasicRandGenerator::C = 0;

BasicRandGenerator::BasicRandGenerator()
{
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
