#ifndef BASICRANDGENERATOR_H
#define BASICRANDGENERATOR_H

#include "randlib_global.h"

/**
 * @brief The BasicRandGenerator class
 * Generate random variable uniformly distributed between 0 and 2^32-1
 */
class RANDLIBSHARED_EXPORT BasicRandGenerator
{
    /// Variables for pseudo generator
    static unsigned long startPoint;
    unsigned long X, Y, Z, W;
    bool C;

public:
    BasicRandGenerator();

    /**
     * @brief getRand
     * Random variable generator Keep It Simply Stupid
     * @return random variable as standard rand() function
     */
    unsigned long getRand();

    /**
     * @brief maxValueInv
     * @return 2^(-32)
     */
    static constexpr double maxInv() { return .23283064e-9; }
};


#endif // BASICRANDGENERATOR_H
