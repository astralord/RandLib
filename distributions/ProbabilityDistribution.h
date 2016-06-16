#ifndef PROBABILITY_DISTRIBUTION_H
#define PROBABILITY_DISTRIBUTION_H

#include <cmath>
#include <string>
#include <utility>

#include "math/RandMath.h"
#include "randlib_global.h"

typedef std::pair <double, double> DoublePair;
typedef std::pair <int, int> IntPair;
typedef std::pair <double, int> DoubleIntPair;

/**
 * @brief The ProbabilityDistribution class
 */
template < typename T >
class RANDLIBSHARED_EXPORT ProbabilityDistribution
{

protected:
    std::string toStringWithPrecision(const double a_value, const int n = 6) const;

public:
    ProbabilityDistribution();
    virtual ~ProbabilityDistribution() {}

    /**
     * @brief name
     * @return title of distribution, for instance "Normal(0, 1)"
     */
    virtual std::string name() const = 0;

    /**
     * @brief variate()
     * @return random variable
     */
    virtual T variate() const = 0;

    /**
     * @brief sample
     * @param outputData
     */
    virtual void sample(std::vector<T> &outputData) const;

    /**
     * @brief F
     * @param x
     * @return P(X < x)
     */
    virtual double F(T x) const = 0;

    /**
     * @brief CumulativeDistributionFunction
     * @param x input vector
     * @param y output vector: y = P(X < x)
     */
    void CumulativeDistributionFunction(const std::vector<T> &x, std::vector<double> &y) const;
};

#endif // PROBABILITY_DISTRIBUTION_H
