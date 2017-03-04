#ifndef PROBABILITY_DISTRIBUTION_H
#define PROBABILITY_DISTRIBUTION_H

#include <string>

#include "math/RandMath.h"
#include "randlib_global.h"

/**
 * @brief The ProbabilityDistribution class
 */
template < typename T >
class RANDLIBSHARED_EXPORT ProbabilityDistribution
{
protected:
    /**
     * @brief MAX_ITER_REJECTION
     * upper boundary for maximum amount of iterations in rejection methods
     * one hundred should be enough to be sure there is a bug
     * (or rejection method is too slow to be used)
     */
    static constexpr double MAX_ITER_REJECTION = 100;

    std::string toStringWithPrecision(const double a_value, const int n = 6) const;

public:
    ProbabilityDistribution();
    virtual ~ProbabilityDistribution() {}

    /**
     * @brief Name
     * @return title of distribution, for instance "Normal(0, 1)"
     */
    virtual std::string Name() const = 0;

    /**
     * @brief MinValue
     * @return minimum possible value that can be achieved by random variable
     */
    virtual T MinValue() const = 0;

    /**
     * @brief MaxValue
     * @return maximum possible value that can be achieved by random variable
     */
    virtual T MaxValue() const = 0;

    /**
     * @brief Variate()
     * @return random variable
     */
    virtual T Variate() const = 0;

    /**
     * @brief Sample
     * @param outputData
     */
    virtual void Sample(std::vector<T> &outputData) const;

    /**
     * @brief F
     * @param x
     * @return P(X < x)
     */
    virtual double F(const T & x) const = 0;

    /**
     * @brief CumulativeDistributionFunction
     * @param x input vector
     * @param y output vector: y = P(X < x)
     */
    void CumulativeDistributionFunction(const std::vector<T> &x, std::vector<double> &y) const;

    /**
     * @brief S
     * @param x
     * @return 1 - P(X < x)
     */
    virtual double S(const T & x) const;

    /**
     * @brief SurvivalFunction
     * @param x input vector
     * @param y output vector: y = 1 - P(X < x)
     */
    void SurvivalFunction(const std::vector<T> &x, std::vector<double> &y) const;

protected:
    enum FIT_ERROR_TYPE {
        WRONG_SAMPLE,
        NOT_APPLICABLE,
        WRONG_RETURN,
        TOO_FEW_ELEMENTS,
        WRONG_LEVEL,
        UNDEFINED_ERROR
    };

    static constexpr char POSITIVITY_VIOLATION[] = "All elements should be positive";
    static constexpr char NON_NEGATIVITY_VIOLATION[] = "All elements should be non-negative";
    static constexpr char UPPER_LIMIT_VIOLATION[] = "No element shouldn be bigger than ";
    static constexpr char LOWER_LIMIT_VIOLATION[] = "No element should be less than ";

    std::string fitError(FIT_ERROR_TYPE fet, const std::string &explanation);
};

#endif // PROBABILITY_DISTRIBUTION_H
