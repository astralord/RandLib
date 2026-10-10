#include "BasicRandGenerator.h"
#include <atomic>
#include <chrono>
#include <exception>
#include <functional>
#include <random>
#include <thread>

unsigned long RandEngine::mix(unsigned long a, unsigned long b, unsigned long c)
{
    a = a - b;  a = a - c;  a = a ^ (c >> 13);
    b = b - c;  b = b - a;  b = b ^ (a << 8);
    c = c - a;  c = c - b;  c = c ^ (b >> 13);
    a = a - b;  a = a - c;  a = a ^ (c >> 12);
    b = b - c;  b = b - a;  b = b ^ (a << 16);
    c = c - a;  c = c - b;  c = c ^ (b >> 5);
    a = a - b;  a = a - c;  a = a ^ (c >> 3);
    b = b - c;  b = b - a;  b = b ^ (a << 10);
    c = c - a;  c = c - b;  c = c ^ (b >> 15);
    return c;
}

unsigned long long RandEngine::splitMix64(unsigned long long &state)
{
    unsigned long long z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

unsigned long RandEngine::getRandomSeed()
{
    static std::atomic<unsigned long long> counter{0};
    unsigned long long entropy = 0;
    try {
        std::random_device device;
        entropy = (static_cast<unsigned long long>(device()) << 32) ^ device();
    } catch (const std::exception &) {
        /// std::random_device may be unavailable; the other sources still differ per call
    }
    entropy ^= static_cast<unsigned long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    entropy ^= static_cast<unsigned long long>(std::hash<std::thread::id>()(std::this_thread::get_id())) << 1;
    entropy += 0x9E3779B97F4A7C15ULL * ++counter;
    return static_cast<unsigned long>(splitMix64(entropy));
}

void JKissRandEngine::Reseed(unsigned long seed)
{
    unsigned long long s = seed;
    X = static_cast<unsigned int>(splitMix64(s));
    do {
        Y = static_cast<unsigned int>(splitMix64(s)); /// xorshift state must not be zero
    } while (Y == 0);
    Z = static_cast<unsigned int>(splitMix64(s));
    C = 1 + static_cast<unsigned int>(splitMix64(s) % 698769067ULL); /// carry in [1, a - 2]
}

unsigned long long JKissRandEngine::Next()
{
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

void JLKiss64RandEngine::Reseed(unsigned long seed)
{
    unsigned long long s = seed;
    X = splitMix64(s);
    do {
        Y = splitMix64(s); /// xorshift state must not be zero
    } while (Y == 0);
    Z1 = static_cast<unsigned int>(splitMix64(s));
    Z2 = static_cast<unsigned int>(splitMix64(s));
    C1 = 1 + static_cast<unsigned int>(splitMix64(s) % 4294584391ULL); /// carries in [1, a - 2]
    C2 = 1 + static_cast<unsigned int>(splitMix64(s) % 4246477507ULL);
}

unsigned long long JLKiss64RandEngine::Next()
{
    X = 1490024343005336237ULL * X + 123456789;
    Y ^= Y << 21;
    Y ^= Y >> 17;
    Y ^= Y << 30;

    unsigned long long t = 4294584393ULL * Z1 + C1;
    C1 = t >> 32;
    Z1 = t;
    t = 4246477509ULL * Z2 + C2;
    C2 = t >> 32;
    Z2 = t;
    return X + Y + Z1 + (static_cast<unsigned long long>(Z2) << 32);
}

void PCGRandEngine::Reseed(unsigned long seed)
{
    Reseed(seed, DEFAULT_SEQUENCE);
}

void PCGRandEngine::Reseed(unsigned long long seed, unsigned long long sequence)
{
    state = 0;
    inc = (sequence << 1u) | 1u;
    Next();
    state += seed;
    Next();
}

unsigned long long PCGRandEngine::Next()
{
    unsigned long long oldstate = state;
    state = oldstate * 6364136223846793005ULL + (inc|1);
    unsigned int xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    unsigned int rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
