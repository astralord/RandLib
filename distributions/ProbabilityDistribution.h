#ifndef PROBABILITY_DISTRIBUTION_H
#define PROBABILITY_DISTRIBUTION_H

#include <cmath>
#include <string>

#include "math/RandMath.h"
#include "randlib_global.h"

/**
 * @brief The double2d struct
 */
struct double2d {
    double x;
    double y;
    double2d & operator=(const double2d & other);
};

/**
 * @brief The ProbabilityDistribution class
 */

template < typename T >
class RANDLIBSHARED_EXPORT ProbabilityDistribution
{

protected:
    std::string toStringWithPrecision(const double a_value, const int n = 6);

public:
    ProbabilityDistribution();
    virtual ~ProbabilityDistribution() {}

    /**
     * @brief name
     * @return title of distribution, for instance "Normal(0, 1)"
     */
    virtual std::string name() = 0;

    /**
     * @brief variate()
     * @return random variable
     */
    virtual T variate() const = 0;

    /**
     * @brief sample
     * @param outputData
     */
    virtual void sample(QVector<T> &outputData) const;

    /**
     * @brief F
     * @param x
     * @return P(X < x)
     */
    virtual double F(T x) const = 0;

    /**
     * @brief cdf
     * @param x input vector
     * @param y output vector: y = P(X < x)
     */
    void cdf(const QVector<T> &x, QVector<double> &y) const;

    /**
     * @brief Mean
     * @return Mathematical expectation
     */
    virtual T Mean() const = 0;
};

#endif // PROBABILITY_DISTRIBUTION_H
