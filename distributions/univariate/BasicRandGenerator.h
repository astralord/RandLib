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
     * @fn getRandomSeed
     * @return seed from std::random_device, mixed with a high-resolution clock,
     * the thread id and a per-process counter
     */
    static unsigned long getRandomSeed();
    /**
     * @fn splitMix64
     * Advances the SplitMix64 state and returns the next output; used to expand a seed
     * into the full state of an engine
     * @param state
     * @return next 64-bit output
     */
    static unsigned long long splitMix64(unsigned long long &state);

public:
    RandEngine() {}
    virtual ~RandEngine() {}

    /**
     * @fn DeriveSeed
     * @param seed seed of a distribution object
     * @param stream index of an auxiliary generator inside that object, from 1
     * @return seed for that auxiliary generator; distinct (seed, stream) pairs give unrelated
     * seeds, unlike seed + stream, which collides with the seed of another object
     */
    static unsigned long DeriveSeed(unsigned long seed, unsigned long stream)
    {
        unsigned long long state = static_cast<unsigned long long>(seed) ^ (0xD1B54A32D192ED03ULL * stream);
        return static_cast<unsigned long>(splitMix64(state));
    }
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
    JKissRandEngine() { this->Reseed(getRandomSeed()); }
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
    JLKiss64RandEngine() { this->Reseed(getRandomSeed()); }
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
    /// Stream selector used by Reseed(seed)
    static constexpr unsigned long long DEFAULT_SEQUENCE = 0xda3e39cb94b95bdbULL;

    PCGRandEngine() { this->Reseed(getRandomSeed()); }
    unsigned long long MinValue() const { return 0; }
    unsigned long long MaxValue() const { return 4294967295UL; }
    void Reseed(unsigned long seed);
    /**
     * @fn Reseed
     * Seeds as the reference pcg32_srandom_r(initstate, initseq)
     * @param seed initial state
     * @param sequence stream selector
     */
    void Reseed(unsigned long long seed, unsigned long long sequence);
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

    /**
     * @fn Bits
     * @param count number of random bits, from 1 to 64
     * @return count random bits, taken from the high end of one or more engine outputs
     */
    unsigned long long Bits(int count)
    {
        const int width = static_cast<int>(maxDecimals());
        if (width >= count)
            return engine.Next() >> (width - count);
        unsigned long long bits = 0;
        int filled = 0;
        while (filled < count) {
            bits = (bits << width) | engine.Next();
            filled += width;
        }
        return bits >> (filled - count);
    }
};

/// Select the engine behind RandGenerator with RANDLIB_JLKISS64 or RANDLIB_PCG32
/// (the CMake option RANDLIB_GENERATOR sets them). JLKISS64RAND is accepted as a legacy alias.
#if defined(RANDLIB_JLKISS64) || defined(JLKISS64RAND)
typedef BasicRandGenerator<JLKiss64RandEngine> RandGenerator;
#elif defined(RANDLIB_PCG32)
typedef BasicRandGenerator<PCGRandEngine> RandGenerator;
#else
typedef BasicRandGenerator<JKissRandEngine> RandGenerator;
#endif


#endif // BASICRANDGENERATOR_H
