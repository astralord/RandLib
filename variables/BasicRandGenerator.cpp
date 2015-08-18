#include "BasicRandGenerator.h"
#include <time.h>

BasicRandGenerator::BasicRandGenerator()
{
}

unsigned long long BasicRandGenerator::variate()
{
    return BasicRandGenerator::JLKISS64();
}

unsigned long long BasicRandGenerator::JLKISS64() {
    static unsigned long long X = 123456789 ^ time(0);
    static unsigned long long Y = X ^ X;
    static unsigned int Z1 = X ^ Y;
    static unsigned int Z2 = X ^ Z1;
    static unsigned int C1 = X ^ Z2;
    static unsigned int C2 = X ^ C1;

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

unsigned long BasicRandGenerator::JKISS32()
{
    static unsigned int X = 123456789 ^ time(0);
    static unsigned int Y = X ^ X;
    static unsigned int Z = X ^ Y;
    static unsigned int W = X ^ Z;
    static bool C = 0;

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
