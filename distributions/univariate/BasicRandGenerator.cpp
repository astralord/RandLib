#include "BasicRandGenerator.h"
#include <time.h>

template <>
unsigned long long BasicRandGenerator<JLKISS64>::variate()
{
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

template <>
unsigned long long BasicRandGenerator<JKISS>::variate()
{
    static unsigned int X = 123456789 ^ time(0);
    static unsigned int C = 6543217;
    static unsigned int Y = X ^ C;
    static unsigned int Z = X ^ Y;
    unsigned long long t = 698769069ULL * Z + C;

    X *= 69069;
    X += 12345;

    Y ^= Y << 13;
    Y ^= Y >> 17;
    Y ^= Y << 5;

    C = t >> 32;
    Z = t;

    return X + Y + Z;
}

template < char Generator >
size_t BasicRandGenerator<Generator>::maxDecimals()
{
    size_t num = 0;
    unsigned long long maxRand = maxValue();
    while (maxRand != 0)
    {
        ++num;
        maxRand >>= 1;
    }
    return num;
}

template class BasicRandGenerator<JLKISS64>;
template class BasicRandGenerator<JKISS>;
