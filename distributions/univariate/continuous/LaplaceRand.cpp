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
        throw std::invalid_argument("Laplace distribution: scale should be positive, but it's equal to "
                                    + std::to_string(scale));
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
long double AsymmetricLaplaceDistribution<RealType>::Entropy() const
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

    if (!RandMath::findRootNewtonFirstOrder<double>([this, sample] (double med)
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
            deviation += this->kappa * (x - this->m);
        else
            deviation -= (x - this->m) / this->kappa;
    }
    deviation /= sample.size();
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
DoublePair AsymmetricLaplaceRand<RealType>::getOneSidedSums(const std::vector<RealType> &sample)
{
    double xPlus = 0.0, xMinus = 0.0;
    for (const RealType & x : sample) {
        if (x < this->m)
            xMinus += (x - this->m);
        else
            xPlus += (x - this->m);
    }
    return DoublePair(xPlus, xMinus);
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::FitAsymmetry(const std::vector<RealType> &sample)
{
    auto [xPlus, xMinus] = getOneSidedSums(sample);
    double gammaN = this->gamma * sample.size();
    double asymmetry = 0;
    if (xPlus == xMinus)
        asymmetry = 1.0;
    else if (xPlus == 0)
        asymmetry = -xMinus / gammaN;
    else if (xMinus == 0) {
        asymmetry = gammaN / xPlus;
    }
    else {
        /// write down coefficients of quartic equation
        double a = xPlus / gammaN;
        double c = (xPlus + xMinus) / gammaN;
        double e = xMinus / gammaN;
        /// find coefficients for solution
        double ae = a * e;
        double delta1 = 2 * std::pow(c, 3) + 9 * c + 27 * (a + e) - 72 * ae * c;
        double delta0 = c * c + 3 + 12 * ae;
        double p2 = delta1 + std::sqrt(delta1 * delta1 - 4 * std::pow(delta0, 3));
        double Q = std::cbrt(0.5 * p2);
        double temp = 8 * a * a;
        double p = (8 * a * c - 3.0) / temp;
        double q = (-1.0 + 4 * a * c + temp) / (a * temp);
        double S = 0.5 * std::sqrt(-2.0 / 3 * p + (Q + delta0 / Q) / (3 * a));
        /// solve by general formula
        double c1 = -0.25 / a, b1 = -4 * S * S - 2 * p, b2 = q / S;
        if (b1 + b2 > 0)
            asymmetry = c1 + S + 0.5 * std::sqrt(b1 + b2);
        else
            asymmetry = c1 - S + 0.5 * std::sqrt(b1 - b2);
    }

    SetAsymmetry(asymmetry);
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
