#include "SumRand.h"

template <typename T1, typename T2>
SumRand<T1, T2>::SumRand(const UnivariateProbabilityDistribution<T1> &leftRV, const UnivariateProbabilityDistribution<T2> &rightRV)
    : X(leftRV), Y(rightRV)
{

}

template <typename T1, typename T2>
double SumRand<T1, T2>::Convolution(const std::function<double (T1)> &funPtrX, const std::function<double (T2)> &funPtrY, T1 x, bool isCDF) const
{
    SUPPORT_TYPE suppX = X.supportType(), suppY = Y.supportType();

    if (suppX != INFINITE_T && suppY != INFINITE_T) {

        bool integrandIsLeftBounded = false, integrandIsRightBounded = false;
        T1 minPoint = 0, maxPoint = 0;

        if (Y.isLeftBounded()) {
            integrandIsLeftBounded = true;
            minPoint = Y.MinValue();
        }

        if (X.isRightBounded() && !isCDF) {
            T1 minPointX = x - X.MaxValue();
            if (integrandIsLeftBounded)
                minPoint = std::max(minPoint, minPointX);
            else {
                minPoint = minPointX;
                integrandIsLeftBounded = true;
            }
        }

        if (Y.isRightBounded()) {
            integrandIsRightBounded = true;
            maxPoint = Y.MaxValue();
        }

        if (X.isLeftBounded()) {
            T1 minPointX = x - X.MinValue();
            if (integrandIsRightBounded)
                maxPoint = std::min(maxPoint, minPointX);
            else {
                maxPoint = minPointX;
                integrandIsRightBounded = true;
            }
        }

        if (minPoint > maxPoint) {
            return 0.0;
        }

        if (integrandIsLeftBounded && integrandIsRightBounded) {
            return Y.ExpectedValue([this, x, funPtrX](double t)
            {
                return funPtrX(x - t);
            }, minPoint, maxPoint);
        }
    }

    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = 0.0;

    if (suppX == INFINITE_T) {
        return Y.ExpectedValue([this, x, funPtrX](double t)
        {
            return funPtrX(x - t);
        }, startPoint);
    }

    return X.ExpectedValue([this, x, funPtrY](double t)
    {
        return funPtrY(x - t);
    }, startPoint);
}


template <typename T1, typename T2>
std::string SumRand<T1, T2>::name() const
{
    return X.name() + " + " + Y.name();
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
    if (this->isLeftBounded() && x < MinValue())
        return 0.0;
    if (this->isRightBounded() && x > MaxValue())
        return 1.0;
    return Convolution([this](T1 t)
    {
        return X.F(t);
    },
    [this](double t)
    {
        return Y.F(t);
    }, x, true);
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
    if (this->isLeftBounded() && x < MinValue())
        return 0.0;
    if (this->isRightBounded() && x > MaxValue())
        return 0.0;
    return Convolution([this](double t)
    {
        return X.f(t);
    },
    [this](double t)
    {
        return Y.f(t);
    }, x, false);
}


double SumDiscreteRand::P(int k) const
{
    if (this->isLeftBounded() && k < MinValue())
        return 0.0;
    if (this->isRightBounded() && k > MaxValue())
        return 0.0;
    return Convolution([this](int x)
    {
        return X.P(x);
    },
    [this](int x)
    {
        return Y.P(x);
    }, k, false);
}

int SumDiscreteRand::Mode() const
{
    // TODO: verify this!!!
    return X.Mode() + Y.Mode();
}

