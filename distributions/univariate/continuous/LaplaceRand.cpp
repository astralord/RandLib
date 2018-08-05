#include "LaplaceRand.h"
#include "../discrete/BernoulliRand.h"

template < typename RealType >
AsymmetricLaplaceDistribution<RealType>::AsymmetricLaplaceDistribution(double shift, double scale, double asymmetry)
    : ShiftedGeometricStableDistribution<RealType>(2.0, 0.0, scale, 0.0, shift)
{
    ShiftedGeometricStableDistribution<RealType>::SetAsymmetry(asymmetry);
    ChangeLocation();
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::ChangeLocation()
{
    this->SetLocation((1.0 - this->kappaSq) * this->gamma * this->kappaInv);
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Laplace distribution: scale should be positive");
    ShiftedGeometricStableDistribution<RealType>::SetScale(scale);
    ChangeLocation();
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::f(const RealType & x) const
{
    return this->pdfLaplace(x - this->m);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::logf(const RealType & x) const
{
    return this->logpdfLaplace(x - this->m);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::F(const RealType & x) const
{
    return this->cdfLaplace(x - this->m);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::S(const RealType & x) const
{
    return this->cdfLaplaceCompl(x - this->m);
}

template < typename RealType >
RealType AsymmetricLaplaceDistribution<RealType>::Variate() const
{
    RealType X = (this->kappa == 1) ? LaplaceRand<RealType>::StandardVariate(this->localRandGenerator)
                                    : AsymmetricLaplaceRand<RealType>::StandardVariate(this->kappa, this->localRandGenerator);
    return this->m + this->gamma * X;
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::Sample(std::vector<RealType> &outputData) const
{
    if (this->kappa == 1) {
        for (RealType & var : outputData)
            var = this->m + this->gamma * LaplaceRand<RealType>::StandardVariate(this->localRandGenerator);
    }
    else {
        for (RealType & var : outputData)
            var = this->m + this->gamma * AsymmetricLaplaceRand<RealType>::StandardVariate(this->kappa, this->localRandGenerator);
    }
}

template < typename RealType >
std::complex<double> AsymmetricLaplaceDistribution<RealType>::CFImpl(double t) const
{
    double bt = this->gamma * t;
    double btSq = bt * bt;
    double denominator = (1 + this->kappaSq * btSq) * (1 + btSq / this->kappaSq);
    std::complex<double> y(std::cos(this->m * t), std::sin(this->m * t));
    std::complex<double> x(1, -this->kappa * bt), z(1, bt * this->kappaInv);
    return x * y * z / denominator;
}

template < typename RealType >
RealType AsymmetricLaplaceDistribution<RealType>::quantileImpl(double p) const
{
    return this->quantileLaplace(p);
}

template < typename RealType >
RealType AsymmetricLaplaceDistribution<RealType>::quantileImpl1m(double p) const
{
    return this->quantileLaplace1m(p);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::Entropy() const
{
    double y = this->kappaInv + this->kappa;
    return std::log1pl(this->gamma * y);
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::FitShift(const std::vector<RealType> &sample)
{
    /// Calculate median (considering asymmetry)
    /// we use root-finding algorithm for median search
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());
    double median = this->GetSampleMean(sample);

    if (!RandMath::findRoot<double>([this, sample] (double med)
    {
        double y = 0.0;
        for (const RealType & x : sample) {
            if (x > med)
                y -= this->kappaSq;
            else if (x < med)
                ++y;
        }
        return y;
    },
    minVar, maxVar, median
    ))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetShift(median);
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::FitScale(const std::vector<RealType> &sample)
{
    double deviation = 0.0;
    for (const RealType & x : sample) {
        if (x > this->m)
            deviation += this->kappaSq * (x - this->m);
        else
            deviation -= (x - this->m);
    }
    deviation /= (this->kappa * sample.size());

    SetScale(deviation);
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::FitShiftAndScale(const std::vector<RealType> &sample)
{
    FitShift(sample);
    FitScale(sample);
}

template class AsymmetricLaplaceDistribution<float>;
template class AsymmetricLaplaceDistribution<double>;
template class AsymmetricLaplaceDistribution<long double>;

template < typename RealType >
String AsymmetricLaplaceRand<RealType>::Name() const
{
    return "Asymmetric-Laplace(" + this->toStringWithPrecision(this->GetShift()) + ", "
                                 + this->toStringWithPrecision(this->GetScale()) + ", "
                                 + this->toStringWithPrecision(this->GetAsymmetry()) + ")";
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::SetAsymmetry(double asymmetry)
{
    ShiftedGeometricStableDistribution<RealType>::SetAsymmetry(asymmetry);
    this->ChangeLocation();
}

template < typename RealType >
RealType AsymmetricLaplaceRand<RealType>::StandardVariate(double asymmetry, RandGenerator &randGenerator)
{
    RealType x = ExponentialRand<RealType>::StandardVariate(randGenerator) / asymmetry;
    RealType y = ExponentialRand<RealType>::StandardVariate(randGenerator) * asymmetry;
    return x - y;
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::FitAsymmetry(const std::vector<RealType> &sample)
{
    double xPlus = 0.0, xMinus = 0.0;
    for (const RealType & x : sample) {
        if (x < this->m)
            xMinus -= (x - this->m);
        else
            xPlus += (x - this->m);
    }

    if (xPlus == xMinus) {
        SetAsymmetry(1.0);
        return;
    }

    double gammaN = this->gamma * sample.size();
    double root = 1.0;
    double minBound, maxBound;
    if (xPlus < -xMinus) { //TODO: why there is minus?
        minBound = 1.0;
        maxBound = std::sqrt(xMinus / xPlus);
    }
    else {
        minBound = std::sqrt(xMinus / xPlus);
        maxBound = 1.0;
    }

    if (!RandMath::findRoot<double>([sample, xPlus, xMinus, gammaN] (double t)
    {
        double tSq = t * t;
        double y = 1.0 - tSq;
        y /= (t * tSq + t);
        y *= gammaN;
        y += xMinus / tSq - xPlus;
        return y;
    }, minBound, maxBound, root))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetAsymmetry(root);
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::FitShiftAndAsymmetry(const std::vector<RealType> &sample)
{
    this->FitShift(sample);
    FitAsymmetry(sample);
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::FitScaleAndAsymmetry(const std::vector<RealType> &sample)
{
    int n = sample.size();
    double xPlus = 0.0, xMinus = 0.0;
    for (const RealType & x : sample) {
        if (x < this->m)
            xMinus -= (x - this->m);
        else
            xPlus += (x - this->m);
    }
    xPlus /= n;
    xMinus /= n;

    if (xMinus == 0) {
        /// X ~ Exp(1 / xPlus) + m
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Sample might be from shifted exponential distribution (no values smaller than shift m)"));
    }
    if (xPlus == 0) {
        /// -X ~ Exp(1 / xMinus) + m
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Sample might be from shifted negative exponential distribution (no values larger than shift m)"));
    }

    double xPlusSqrt = std::sqrt(xPlus), xMinusSqrt = std::sqrt(xMinus);
    double scale = xPlusSqrt + xMinusSqrt;
    scale *= std::sqrt(xPlusSqrt * xMinusSqrt);

    this->SetScale(scale);
    SetAsymmetry(std::pow(xMinus / xPlus, 0.25));
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::Fit(const std::vector<RealType> &sample)
{
    this->FitShift(sample);
    FitScaleAndAsymmetry(sample);
}

template class AsymmetricLaplaceRand<float>;
template class AsymmetricLaplaceRand<double>;
template class AsymmetricLaplaceRand<long double>;

template < typename RealType >
String LaplaceRand<RealType>::Name() const
{
    return "Laplace(" + this->toStringWithPrecision(this->GetShift()) + ", "
                      + this->toStringWithPrecision(this->GetScale()) + ")";
}

template < typename RealType >
RealType LaplaceRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    RealType W = ExponentialRand<RealType>::StandardVariate(randGenerator);
    return BernoulliRand::StandardVariate(randGenerator) ? W : -W;
}


template class LaplaceRand<float>;
template class LaplaceRand<double>;
template class LaplaceRand<long double>;
