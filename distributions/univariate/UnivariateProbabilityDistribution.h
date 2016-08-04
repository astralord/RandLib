#ifndef UNIVARIATEPROBABILITYDISTRIBUTION_H
#define UNIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"

enum SUPPORT_TYPE {
    FINITE_T,
    RIGHTSEMIFINITE_T,
    LEFTSEMIFINITE_T,
    INFINITE_T
};

/**
 * @brief The UnivariateProbabilityDistribution class
 */
template < typename T >
class RANDLIBSHARED_EXPORT UnivariateProbabilityDistribution : public ProbabilityDistribution<T>
{
public:
    UnivariateProbabilityDistribution();
    virtual ~UnivariateProbabilityDistribution() {}

    /**
     * @brief supportType
     * @return type of support
     */
    virtual SUPPORT_TYPE supportType() const = 0;

    /**
     * @brief isLeftBounded
     * @return true if distribution is bounded from the left
     */
    bool isLeftBounded() const;

    /**
     * @brief isRightBounded
     * @return true if distribution is bounded from the right
     */
    bool isRightBounded() const;

    /**
     * @brief MinValue
     * @return minimum possible value that can be achieved by random variable
     */
    virtual T MinValue() const override = 0;

    /**
     * @brief MaxValue
     * @return maximum possible value that can be achieved by random variable
     */
    virtual T MaxValue() const override = 0;

    /**
     * @brief Mean
     * @return Mathematical expectation
     */
    virtual double Mean() const = 0;

    /**
     * @brief Variance
     * @return Variance of random variable
     */
    virtual double Variance() const = 0;

    /**
     * @brief Quantile
     * @param p
     * @return such x that F(x) = p
     */
    virtual double Quantile(double p) const = 0;

    /**
     * @brief QuantileFunction
     * @param p
     * @return fills vector y with Quantile(p)
     */
    void QuantileFunction(const std::vector<double> &p, std::vector<double> &y);

    /**
     * @brief CF
     * @param x
     * @return Characteristic function (inverse Fourier transform of probability function)
     */
    virtual std::complex<double> CF(double t) const;

    /**
     * @brief CharacteristicFunction
     * @param x input vector
     * @param y output vector: y = CF(x)
     */
    void CharacteristicFunction(const std::vector<double> &t, std::vector<std::complex<double>> &y) const;

    /**
     * @brief Hazard
     * return hazard function: pdf (or pmf) / (1 - cdf)
     * @param x input parameter
     */
    virtual double Hazard(double x) const = 0;

    /**
     * @brief HazardFunction
     * @param x input vector
     * @param y output vector: y = Hazard(x)
     */
    void HazardFunction(const std::vector<double> &x, std::vector<double> &y) const;

    /**
     * @brief ExpectedValue
     * @param funPtr pointer on function g(x) which expected value should be returned
     * @param startPoint argument in which vicinity value of g(x) definitely wouldn't be zero
     * @return E[g(x)]
     */
    virtual double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const = 0;

    /**
     * @brief ExpectedValue
     * @param funPtr pointer on function g(x) with finite support which expected value should be returned
     * @param minPoint min{x | g(x) != 0}
     * @param maxPoint max{x | g(x) != 0}
     * @return E[g(x)]
     */
    virtual double ExpectedValue(const std::function<double (double)> &funPtr, T minPoint, T maxPoint) const = 0;

    /**
     * @brief Median
     * @return such x that F(x) = 0.5
     */
    virtual double Median() const;

    /**
     * @brief Mode
     * @return the most probable value
     */
    virtual T Mode() const = 0;

    /**
     * @brief Skewness
     * @return E[((X - μ) / σ) ^ 3]
     * where mu is central moment and sigma is standard deviation
     */
    virtual double Skewness() const;

    /**
     * @brief Kurtosis
     * @return unbiased kurtosis = μ_4 / σ ^ 4
     */
    virtual double Kurtosis() const;

    /**
     * @brief ExcessKurtosis
     * @return E[((X - μ) / σ) ^ 4]  - 3
     * (fourth moment around the mean divided by the square of the variance of the probability distribution minus 3)
     */
    virtual double ExcessKurtosis() const;

    /**
     * @brief SecondMoment
     * @return E[X^2]
     */
    virtual double SecondMoment() const;

    /**
     * @brief ThirdMoment
     * @return E[X^3]
     */
    virtual double ThirdMoment() const;

    /**
     * @brief FourthMoment
     * @return E[X^4]
     */
    virtual double FourthMoment() const;

    /**
     * @brief Likelihood
     * @param sample
     * @return product of f(x_i)
     */
    virtual double Likelihood(const std::vector<T> &sample) const = 0;

    /**
     * @brief LogLikelihood
     * @param sample
     * @return sum of log(f(x_i))
     */
    virtual double LogLikelihood(const std::vector<T> &sample) const = 0;
};

#endif // UNIVARIATEPROBABILITYDISTRIBUTION_H
