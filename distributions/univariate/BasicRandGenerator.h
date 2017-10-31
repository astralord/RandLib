#ifndef BASICRANDGENERATOR_H
#define BASICRANDGENERATOR_H

#include "RandLib_global.h"
#include <time.h>

enum GeneratorType {
    JKISS, ///< period is 2^1271
    JLKISS64 ///< period is 2^250
};


/**
 * @brief The BasicRandGenerator class
 * Class for generators of random number, spreaded uniformly
 */
template < char Generator >
class RANDLIBSHARED_EXPORT BasicRandGenerator
{

public:
    BasicRandGenerator() {}

    static unsigned long long Variate();

    static constexpr unsigned long long MinValue() {
        return 0;
    }

    static constexpr unsigned long long MaxValue() {
        return (Generator == JLKISS64) ? 18446744073709551615ULL : 4294967295UL;
    }
    static size_t maxDecimals();
};

#ifdef JLKISS64RAND
typedef BasicRandGenerator<JLKISS64> RandGenerator;
#else
typedef BasicRandGenerator<JKISS> RandGenerator;
#endif


#endif // BASICRANDGENERATOR_H
