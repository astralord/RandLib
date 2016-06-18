#include "SumRand.h"

template <typename T1, typename T2>
SumRand<T1, T2>::SumRand(const UnivariateProbabilityDistribution<T1> &leftRV, const UnivariateProbabilityDistribution<T2> &rightRV)
    : X(leftRV), Y(rightRV)
{

}

template <typename T1, typename T2>
std::string SumRand<T1, T2>::name() const
{
    return X.name() + "+" + Y.name();
}

template <typename T1, typename T2>
SUPPORT_TYPE SumRand<T1, T2>::supportType() const
{
    SUPPORT_TYPE XT = X.supportType(), YT = Y.supportType();
    if (XT == FINITE_T)
        return YT;
    if (YT == FINITE_T)
        return XT;
    if (XT == INFINITE_T || YT == INFINITE_T)
        return INFINITE_T;
    if (XT == RIGHTSEMIFINITE_T)
        return (YT == RIGHTSEMIFINITE_T) ? RIGHTSEMIFINITE_T : INFINITE_T;
    return (YT == RIGHTSEMIFINITE_T) ? INFINITE_T : LEFTSEMIFINITE_T;
}

template <typename T1, typename T2>
double SumRand<T1, T2>::MinValue() const
{
    return X.MinValue() + Y.MinValue();
}

template <typename T1, typename T2>
double SumRand<T1, T2>::MaxValue() const
{
    return X.MaxValue() + Y.MaxValue();
}

template <typename T1, typename T2>
double SumRand<T1, T2>::F(T1 x) const
{
    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = 0.0;

    if (X.supportType() == FINITE_T &&
            (Y.supportType() != FINITE_T || X.MaxValue() - X.MinValue() <= Y.MaxValue() - Y.MinValue())) {
        return X.ExpectedValue([this, x](double t)
        {
            return Y.F(x - t);
        }, startPoint);
    }

    return Y.ExpectedValue([this, x](double t)
    {
        return X.F(x - t);
    }, startPoint);
}

template <typename T1, typename T2>
T1 SumRand<T1, T2>::variate() const
{
    return X.variate() + Y.variate();
}

template <typename T1, typename T2>
double SumRand<T1, T2>::Mean() const
{
    return X.Mean() + Y.Mean();
}

template <typename T1, typename T2>
double SumRand<T1, T2>::Variance() const
{
    return X.Variance() + Y.Variance();
}

template <typename T1, typename T2>
std::complex<double> SumRand<T1, T2>::CF(double t) const
{
    return X.CF(t) * Y.CF(t);
}

template <typename T1, typename T2>
double SumRand<T1, T2>::Skewness() const
{
    double skewX = X.Skewness();
    double skewY = Y.Skewness();
    double varX = X.Variance();
    double varY = Y.Variance();
    double skew = skewX * std::pow(varX, 1.5);
    skew += skewY * std::pow(varY, 1.5);
    return skew / std::pow(varX + varY, 1.5);
}

template <typename T1, typename T2>
double SumRand<T1, T2>::ExcessKurtosis() const
{
    double varXSq = X.Variance();
    varXSq *= varXSq;
    double varYSq = Y.Variance();
    varYSq *= varYSq;
    double kurtX = X.Kurtosis();
    double kurtY = Y.Kurtosis();
    double kurt = kurtX * varXSq * varXSq;
    kurt += kurtY * varYSq * varYSq;
    kurt += 6 * varXSq * varYSq;
    double denom = varXSq + varYSq;
    return kurt / (denom * denom);
}

double SumContinuousRand::f(double x) const
{
    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = 0.0;
    if (X.supportType() == FINITE_T &&
            (Y.supportType() != FINITE_T || X.MaxValue() - X.MinValue() <= Y.MaxValue() - Y.MinValue())) {
        return X.ExpectedValue([this, x](double t)
        {
            return Y.f(x - t);
        }, startPoint);
    }

    return Y.ExpectedValue([this, x](double t)
    {
        return X.f(x - t);
    }, startPoint);
}


double SumDiscreteRand::P(int k) const
{
    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = 0.0;
    if (X.supportType() == FINITE_T &&
            (Y.supportType() != FINITE_T || X.MaxValue() - X.MinValue() <= Y.MaxValue() - Y.MinValue())) {
        return X.ExpectedValue([this, k](double t)
        {
            return Y.P(k - t);
        }, startPoint);
    }

    return Y.ExpectedValue([this, k](double t)
    {
        return X.P(k - t);
    }, startPoint);
}

int SumDiscreteRand::Mode() const
{
    // TODO: verify this!!!
    return X.Mode() + Y.Mode();
}
