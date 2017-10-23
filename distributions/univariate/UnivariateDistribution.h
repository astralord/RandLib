#ifndef UNIVARIATEDISTRIBUTION_H
#define UNIVARIATEDISTRIBUTION_H

#include "../ProbabilityDistribution.h"

enum SUPPORT_TYPE {
    FINITE_T,
    RIGHTSEMIFINITE_T,
    LEFTSEMIFINITE_T,
    INFINITE_T
};

/**
 * @brief The UnivariateDistribution class <BR>
 * Abstract class for all univariate probability distributions
 */
template < typename T >
class RANDLIBSHARED_EXPORT UnivariateDistribution : public ProbabilityDistribution<T>
{
protected:
    UnivariateDistribution();
    virtual ~UnivariateDistribution() {}

public:
    /**
     * @fn SupportType
     * @return type of support
     */
    virtual SUPPORT_TYPE SupportType() const = 0;

    /**
     * @fn isLeftBounded
     * @return true if distribution is bounded from the left
     */
    bool isLeftBounded() const;

    /**
     * @fn isRightBounded
     * @return true if distribution is bounded from the right
     */
    bool isRightBounded() const;

    /**
     * @fn Mean
     * @return Mathematical expectation
     */
    virtual double Mean() const = 0;

    /**
     * @fn Variance
     * @return Variance of random variable
     */
    virtual double Variance() const = 0;

private:
    /**
     * @fn quantileImpl
     * @param p
     * @return such x that F(x) = p
     */
    virtual T quantileImpl(double p) const = 0;

    /**
     * @fn quantileImpl1m
     * @param p
     * @return such x that F(x) = 1 - p
     */
    virtual T quantileImpl1m(double p) const = 0;

protected:
    /**
     * @fn CFimpl
     * @param t
     * @return Characteristic function (inverse Fourier transform of pdf)
     */
    virtual std::complex<double> CFImpl(double t) const;

    /**
     * @fn ExpectedValue
     * @param funPtr pointer on function g(x) with finite support which expected value should be returned
     * @param minPoint min{x | g(x) ≠ 0}
     * @param maxPoint max{x | g(x) ≠ 0}
     * @return E[g(x)]
     */
    virtual double ExpectedValue(const std::function<double (double)> &funPtr, T minPoint, T maxPoint) const = 0;
public:
    /**
     * @fn Quantile
     * @param p
     * @return quantileImpl(p) if p is in (0, 1)
     */
    double Quantile(double p) const;

    /**
     * @fn Quantile1m
     * @param p
     * @return quantileImpl1m(p) if p is in (0, 1)
     */
    double Quantile1m(double p) const;

    /**
     * @fn QuantileFunction
     * @param p
     * @return fills vector y with Quantile(p)
     */
    void QuantileFunction(const std::vector<double> &p, std::vector<double> &y);

    /**
     * @fn CF
     * @param t
     * @return CFImpl for t != 0 and 1 otherwise
     */
    std::complex<double> CF(double t) const;

    /**
     * @fn CharacteristicFunction
     * @param x input vector
     * @param y output vector: y = CF(x)
     */
    void CharacteristicFunction(const std::vector<double> &t, std::vector<std::complex<double>> &y) const;

    /**
     * @fn Hazard
     * @param x input parameter
     * @return hazard function
     */
    virtual double Hazard(double x) const = 0;

    /**
     * @fn HazardFunction
     * @param x input vector
     * @param y output vector: y = Hazard(x)
     */
    void HazardFunction(const std::vector<double> &x, std::vector<double> &y) const;

    /**
     * @fn Median
     * @return such x that F(x) = 0.5
     */
    virtual T Median() const;

    /**
     * @fn Mode
     * @return the most probable value
     */
    virtual T Mode() const = 0;

    /**
     * @fn Skewness
     * @return E[((X - μ) / σ) ^ 3]
     * where μ is central moment and σ is standard deviation
     */
    virtual double Skewness() const;

    /**
     * @fn Kurtosis
     * @return unbiased kurtosis = μ_4 / σ ^ 4
     */
    virtual double Kurtosis() const;

    /**
     * @fn ExcessKurtosis
     * @return E[((X - μ) / σ) ^ 4]  - 3
     * (fourth moment around the mean divided by the square of the variance of the probability distribution minus 3)
     */
    virtual double ExcessKurtosis() const;

    /**
     * @fn SecondMoment
     * @return E[X^2]
     */
    virtual double SecondMoment() const;

    /**
     * @fn ThirdMoment
     * @return E[X^3]
     */
    virtual double ThirdMoment() const;

    /**
     * @fn FourthMoment
     * @return E[X^4]
     */
    virtual double FourthMoment() const;

    /**
     * @fn LikelihoodFunction
     * @param sample
     * @return likelihood function for given sample
     */
    virtual double LikelihoodFunction(const std::vector<T> &sample) const = 0;

    /**
     * @fn LogLikelihoodFunction
     * @param sample
     * @return logarithm of likelihood function for given sample
     */
    virtual double LogLikelihoodFunction(const std::vector<T> &sample) const = 0;

protected:
    /**
     * @fn allElementsAreNotBiggerThan
     * @param value
     * @param sample
     * @return true if all elements in sample are not bigger than value
     */
    static bool allElementsAreNotBiggerThan(T value, const std::vector<T> &sample);

    /**
     * @fn allElementsAreNotLessThan
     * @param value
     * @param sample
     * @return true if all elements in sample are not less than value
     */
    static bool allElementsAreNotLessThan(T value, const std::vector<T> &sample);

    /**
     * @fn allElementsAreNonNegative
     * @param sample
     * @return true if all elements in sample are non-negative
     */
    static bool allElementsAreNonNegative(const std::vector<T> &sample);

    /**
     * @fn allElementsArePositive
     * @param sample
     * @return true if all elements in sample are positive
     */
    static bool allElementsArePositive(const std::vector<T> &sample);

public:
    /**
     * @fn GetSampleSum
     * @param sample
     * @return sum of all elements in a sample
     */
    static double GetSampleSum(const std::vector<T> &sample);

    /**
     * @fn GetSampleMean
     * @param sample
     * @return arithmetic average
     */
    static double GetSampleMean(const std::vector<T> &sample);

    /**
     * @fn GetSampleVariance
     * @param sample
     * @param known mean
     * @return sample second central moment
     */
    static double GetSampleVariance(const std::vector<T> &sample, double mean);

    /**
     * @brief GetSampleMeanAndVariance
     * @param sample
     * @return sample mean and variance
     */
    static DoublePair GetSampleMeanAndVariance(const std::vector<T> &sample);

    /**
     * @brief GetSampleStatistics
     * @param sample
     * @return sample mean, variance, skewness and excess kurtosis
     */
    static std::tuple<double, double, double, double> GetSampleStatistics(const std::vector<T> &sample);
};

#endif // UNIVARIATEDISTRIBUTION_H
