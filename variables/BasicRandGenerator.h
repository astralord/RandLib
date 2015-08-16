#ifndef BASICRANDGENERATOR_H
#define BASICRANDGENERATOR_H

#include "randlib_global.h"

enum GeneratorType {
    FASTKISS, // 2 ^ 121
    JLKISS64, // 2 ^ 250
    MWC256, // 2 ^ 8222
    CMWC4096, // 2 ^ 131086
    SUPERKISS // 54767 * 2 ^ 1337279
};

/**
 * @brief The BasicRandGenerator class
 */
class RANDLIBSHARED_EXPORT BasicRandGenerator
{
    /// Variables for pseudo generator
    static unsigned long long startPoint;
    static unsigned long long X, Y;
    static unsigned int Z1, Z2, C1, C2;

public:
    BasicRandGenerator();

    static unsigned long long variate();

    static constexpr unsigned long long maxValue() { return 18446744073709551615ULL; }
};


#endif // BASICRANDGENERATOR_H
