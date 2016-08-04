#ifndef SUMRAND_H
#define SUMRAND_H

#include "UnivariateProbabilityDistribution.h"
#include "continuous/ContinuousDistribution.h"
#include "discrete/DiscreteDistribution.h"

/**
 * @brief The SumRand class
 */
template <typename T1, typename T2>
class RANDLIBSHARED_EXPORT SumRand : public virtual UnivariateProbabilityDistribution<T1>
{
    const UnivariateProbabilityDistribution<T1> &X;
    const UnivariateProbabilityDistribution<T2> &Y;

protected:
    double Convolution(const std::function<double (T1)> &funPtrX, const std::function<double (T2)> &funPtrY, T1 x, bool isCDF) const;

public:
    SumRand(const UnivariateProbabilityDistribution<T1> & leftRV, const UnivariateProbabilityDistribution<T2> & rightRV);
    virtual ~SumRand() {}
    std::string name() const override;

    SUPPORT_TYPE supportType() const override;
    T1 MinValue() const override;
    T1 MaxValue() const override;

    double F(T1 x) const override;
    T1 variate() const override;

    double Mean() const override;
    double Variance() const override;
    std::complex<double> CF(double t) const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The SumContinuousRand class
 */
class RANDLIBSHARED_EXPORT SumContinuousRand : public SumRand<double, double>, public ContinuousDistribution
{
    const ContinuousDistribution &X, &Y;
public:
    SumContinuousRand(const ContinuousDistribution & leftRV, const ContinuousDistribution & rightRV) :
        SumRand<double, double>(leftRV, rightRV),
        X(leftRV), Y(rightRV)
    {}
    virtual ~SumContinuousRand() {}

    double f(double x) const override;

    friend const SumContinuousRand operator+(const ContinuousDistribution& left, const ContinuousDistribution& right);
};

const SumContinuousRand operator+(const ContinuousDistribution& left, const ContinuousDistribution& right) {
    return SumContinuousRand(left, right);
}


/**
 * @brief The SumDiscreteRand class
 */
class RANDLIBSHARED_EXPORT SumDiscreteRand : public SumRand<int, int>, public DiscreteDistribution
{
    const DiscreteDistribution &X, &Y;
public:
    SumDiscreteRand(const DiscreteDistribution & leftRV, const DiscreteDistribution & rightRV) :
        SumRand<int, int>(leftRV, rightRV),
        X(leftRV), Y(rightRV)
    {}
    virtual ~SumDiscreteRand() {}

    double P(int k) const override;

    friend const SumDiscreteRand operator+(const DiscreteDistribution& left, const DiscreteDistribution& right);
};

const SumDiscreteRand operator+(const DiscreteDistribution& left, const DiscreteDistribution& right) {
    return SumDiscreteRand(left, right);
}

/**
 * @brief The SumContinuousDiscreteRand class
 */
class RANDLIBSHARED_EXPORT SumContinuousDiscreteRand : public SumRand<double, int>, public ContinuousDistribution
{
    const ContinuousDistribution &X;
    const DiscreteDistribution &Y;
public:
    SumContinuousDiscreteRand(const ContinuousDistribution & leftRV, const DiscreteDistribution & rightRV) :
        SumRand<double, int>(leftRV, rightRV),
        X(leftRV), Y(rightRV)
    {}
    virtual ~SumContinuousDiscreteRand() {}

    double f(double x) const override;

    friend const SumContinuousDiscreteRand operator+(const ContinuousDistribution& left, const DiscreteDistribution& right);
    friend const SumContinuousDiscreteRand operator+(const DiscreteDistribution& left, const ContinuousDistribution& right);
};

const SumContinuousDiscreteRand operator+(const ContinuousDistribution& left, const DiscreteDistribution& right) {
    return SumContinuousDiscreteRand(left, right);
}

const SumContinuousDiscreteRand operator+(const DiscreteDistribution& left, const ContinuousDistribution& right) {
    return SumContinuousDiscreteRand(right, left);
}

#endif // SUMRAND_H
