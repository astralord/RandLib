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
    virtual ~RandEngine() {}
    virtual unsigned long long MinValue() const = 0;
    virtual unsigned long long MaxValue() const = 0;
    virtual void Reseed(unsigned long seed) = 0;
    virtual unsigned long long Next() = 0;
};

/**
 * @brief The JKissRandEngine class
 */
class RANDLIBSHARED_EXPORT JKissRandEngine : public RandEngine
{
    unsigned int X{};
    unsigned int C{};
    unsigned int Y{};
    unsigned int Z{};

public:
    JKissRandEngine() { this->Reseed(getSeed()); }
    unsigned long long MinValue() const { return 0; }
    unsigned long long MaxValue() const { return 4294967295UL; }
    void Reseed(unsigned long seed);
    unsigned long long Next();
};

/**
 * @brief The JLKiss64RandEngine class
 */
class RANDLIBSHARED_EXPORT JLKiss64RandEngine : public RandEngine
{
    unsigned long long X{};
    unsigned long long Y{};
    unsigned int Z1{};
    unsigned int Z2{};
    unsigned int C1{};
    unsigned int C2{};

public:
    JLKiss64RandEngine() { this->Reseed(getSeed()); }
    unsigned long long MinValue() const { return 0; }
    unsigned long long MaxValue() const { return 18446744073709551615ULL; }
    void Reseed(unsigned long seed);
    unsigned long long Next();
};

/**
 * @brief The PCGRandEngine class
 * Random number generator, taken from http://www.pcg-random.org/
 */
class RANDLIBSHARED_EXPORT PCGRandEngine : public RandEngine
{
    unsigned long long state{};
    unsigned long long inc{};

public:
    PCGRandEngine() { this->Reseed(getSeed()); }
    unsigned long long MinValue() const { return 0; }
    unsigned long long MaxValue() const { return 4294967295UL; }
    void Reseed(unsigned long seed);
    unsigned long long Next();
};

/**
 * @brief The BasicRandGenerator class
 * Class for generators of random number, evenly spreaded from 0 to some integer value
 */
template <class Engine>
class RANDLIBSHARED_EXPORT BasicRandGenerator
{
    static_assert(std::is_base_of<RandEngine, Engine>::value, "Engine must be a descendant of RandEngine");

    Engine engine{};

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

    unsigned long long Variate() { return engine.Next(); }
    size_t maxDecimals() { return getDecimals(engine.MaxValue()); }
    unsigned long long MaxValue() { return engine.MaxValue(); }
    void Reseed(unsigned long seed) { engine.Reseed(seed); }
};

#ifdef JLKISS64RAND
typedef BasicRandGenerator<JLKiss64RandEngine> RandGenerator;
#else
typedef BasicRandGenerator<PCGRandEngine> RandGenerator;
#endif


#endif // BASICRANDGENERATOR_H
