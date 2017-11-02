#ifndef BASICRANDGENERATOR_H
#define BASICRANDGENERATOR_H

#include "RandLib_global.h"
#include <type_traits>
#include <cstddef>

/**
 * @brief The RandEngine class
 */
class RANDLIBSHARED_EXPORT RandEngine
{
protected:
    /**
     * @fn mix
     * Robert Jenkins' 96 bit mix function
     * @param a
     * @param b
     * @param c
     * @return mixing part of the hash function
     */
    static unsigned long mix(unsigned long a, unsigned long b, unsigned long c);
    /**
     * @fn getSeed
     * @return
     */
    static unsigned long getSeed();

public:
    RandEngine() {}
    static constexpr unsigned long long MinValue() { return 0; }
};

/**
 * @brief The JKissRandEngine class
 */
class RANDLIBSHARED_EXPORT JKissRandEngine : public RandEngine
{
    static thread_local unsigned int X;
    static thread_local unsigned int C;
    static thread_local unsigned int Y;
    static thread_local unsigned int Z;

public:
    JKissRandEngine() { Reseed(getSeed()); }
    void Reseed(unsigned long seed);
    static unsigned long long Next();
    static constexpr unsigned long long MaxValue() { return 4294967295UL; }
};

/**
 * @brief The JLKiss64RandEngine class
 */
class RANDLIBSHARED_EXPORT JLKiss64RandEngine : public RandEngine
{
    static thread_local unsigned long long X;
    static thread_local unsigned long long Y;
    static thread_local unsigned int Z1;
    static thread_local unsigned int Z2;
    static thread_local unsigned int C1;
    static thread_local unsigned int C2;

public:
    JLKiss64RandEngine() { Reseed(getSeed()); }
    void Reseed(unsigned long seed);
    static unsigned long long Next();
    static constexpr unsigned long long MaxValue() { return 18446744073709551615ULL; }
};

/**
 * @brief The BasicRandGenerator class
 * Class for generators of random number, evenly spreaded from 0 to some integer value
 */
template < class Engine >
class RANDLIBSHARED_EXPORT BasicRandGenerator
{
    static_assert(std::is_base_of<RandEngine, Engine>::value, "Engine must inherit from RandEngine");

    /**
     * @fn getDecimals
     * @param value
     * @return decimals of given value
     */
    static size_t getDecimals(unsigned long long value)
    {
        size_t num = 0;
        unsigned long long maxRand = value;
        while (maxRand != 0) {
            ++num;
            maxRand >>= 1;
        }
        return num;
    }

public:
    BasicRandGenerator() {}

    static unsigned long long Variate() { return Engine::Next(); }
    static size_t maxDecimals() { return getDecimals(Engine::MaxValue()); }
    static constexpr unsigned long long MaxValue() { return Engine::MaxValue(); }
};

#ifdef JLKISS64RAND
typedef BasicRandGenerator<JLKiss64RandEngine> RandGenerator;
#else
typedef BasicRandGenerator<JKissRandEngine> RandGenerator;
#endif


#endif // BASICRANDGENERATOR_H
