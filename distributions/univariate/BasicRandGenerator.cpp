#include "BasicRandGenerator.h"

template <>
unsigned long long BasicRandGenerator<JLKISS64>::Variate()
{
    static unsigned long long X = 123456789123ULL ^ time(0);
    static unsigned long long Y = 987654321987ULL;
    static unsigned int Z1 = 43219876;
    static unsigned int Z2 = 6543217;
    static unsigned int C1 = 21987643;
    static unsigned int C2 = 1732654;

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
    return X + Y + Z1 + (static_cast<unsigned long long>(Z2) << 32);
}

template <>
unsigned long long BasicRandGenerator<JKISS>::Variate()
{
    static unsigned int X = 123456789 ^ time(0);
    static unsigned int C = 6543217;
    static unsigned int Y = 987654321;
    static unsigned int Z = 43219876;
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
    unsigned long long maxRand = MaxValue();
    while (maxRand != 0)
    {
        ++num;
        maxRand >>= 1;
    }
    return num;
}

template class BasicRandGenerator<JLKISS64>;
template class BasicRandGenerator<JKISS>;
