#ifndef UNIVARIATEPROBABILITYDISTRIBUTION_H
#define UNIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"

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
     * @brief quantile
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
     * @brief ExpectedValue
     * @param funPtr function which expected value should be returned
     * @param startPoint argument in which vicinity value of funPtr definitely wouldn't be zero
     * @return E[funPtr(x)]
     */
    virtual double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const = 0;

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
     * @return E[((X - mu) / sigma) ^ 3]
     * where mu is central moment and sigma is standard deviation
     */
    virtual double Skewness() const;

    /**
     * @brief Kurtosis
     * @return unbiased kurtosis = mu_4 / sigma ^ 4
     */
    virtual double Kurtosis() const;

    /**
     * @brief ExcessKurtosis
     * @return E[((X - mu) / sigma) ^ 4]  - 3
     * (fourth moment around the mean divided by the square of the variance of the probability distribution minus 3)
     */
    virtual double ExcessKurtosis() const;
};

#endif // UNIVARIATEPROBABILITYDISTRIBUTION_H
