#ifndef BASICRANDGENERATOR_H
#define BASICRANDGENERATOR_H

#include "randlib_global.h"

enum GeneratorType {
    JKISS32, // 2 ^ 121
    JLKISS64, // 2 ^ 250
    /*MWC256, // 2 ^ 8222
    CMWC4096, // 2 ^ 131086
    SUPERKISS // 54767 * 2 ^ 1337279*/
};


/**
 * @brief The BasicRandGenerator class
 */
template < char Generator = JLKISS64>
class RANDLIBSHARED_EXPORT BasicRandGenerator
{

public:
    BasicRandGenerator() {}

    static unsigned long long variate();
    static unsigned long long rand_JLKISS64();
    static unsigned long rand_JKISS32();

    static constexpr unsigned long long maxValue() {
        return (Generator == JLKISS64) ? 18446744073709551615ULL : 4294967295UL;
    }
    static size_t maxDecimals();
};

#ifdef JKISS32RAND
typedef BasicRandGenerator<JKISS32> RandGenerator;
#else
typedef BasicRandGenerator<> RandGenerator;
#endif


#endif // BASICRANDGENERATOR_H
