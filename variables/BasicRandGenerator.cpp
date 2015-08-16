#include "BasicRandGenerator.h"
#include <time.h>

unsigned long long BasicRandGenerator::startPoint = 123456789 ^ time(0);
unsigned long long BasicRandGenerator::X = startPoint;
unsigned long long BasicRandGenerator::Y = startPoint ^ X;
unsigned int BasicRandGenerator::Z1 = startPoint ^ Y;
unsigned int BasicRandGenerator::Z2 = startPoint ^ Z1;
unsigned int BasicRandGenerator::C1 = startPoint ^ Z2;
unsigned int BasicRandGenerator::C2 = startPoint ^ C1;

BasicRandGenerator::BasicRandGenerator()
{
}

unsigned long long BasicRandGenerator::variate()
{
    /*
    Y ^= Y << 5;
    Y ^= Y >> 7;
    Y ^= Y << 22;

    int t = Z + W + C;
    Z = W;
    C = t < 0;
    W = t & 2147483647;
    X += 1411392427;

    return X + Y + W;*/


    unsigned long long t;

     X = 1490024343005336237ULL * X + 123456789;
     Y ^= Y << 21;
     Y ^= Y >> 17;
     Y ^= Y << 30;
     t = 4294584393ULL * Z1 + C1;
     C1 = t >> 32;
     Z1 = t;
     t = 4246477509ULL * Z2 + C2;
     C2 = t >> 32;
     Z2 = t;
     return X + Y + Z1 + ((unsigned long long)Z2 << 32);

}
