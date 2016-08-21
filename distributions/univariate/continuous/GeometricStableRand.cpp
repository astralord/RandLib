#include "GeometricStableRand.h"
#include "LaplaceRand.h"
#include "LevyRand.h"
#include "UniformRand.h"
#include "CauchyRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location)
    : LimitingDistribution(exponent, skewness, scale, location),
      Z(exponent, skewness)
{
    SetParameters(exponent, skewness, scale, location);
}

std::string GeometricStableRand::Name() const
{
    return "Geometric Stable("
            + toStringWithPrecision(GetExponent()) + ", "
            + toStringWithPrecision(GetSkewness()) + ", "
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}

void GeometricStableRand::SetParameters(double exponent, double skewness, double scale, double location)
{
    LimitingDistribution::SetParameters(exponent, skewness);
    LimitingDistribution::SetScale(scale);
    LimitingDistribution::SetLocation(location);

    Z.SetParameters(alpha, beta);

    if (alpha == 2) {
        double sigma2 = sigma + sigma;
        k = mu * mu + sigma2 * sigma2;
        k = std::sqrt(k) - mu;
        k /= sigma2;
        kInv = 1.0 / k;
        kSq = k * k;
        pdfCoef = 1.0 / (sigma * (k + kInv));
        cdfCoef = -std::log1p(kSq);
    }
}

double GeometricStableRand::pdfLaplace(double x) const
{
    double y = x / sigma;
    y *= (x < 0) ? kInv : -k;
    y = std::exp(y);
    return pdfCoef * y;
}

double GeometricStableRand::cdfLaplace(double x) const
{
    double y = x / sigma;
    if (x < 0) {
        y *= kInv;
        y = std::exp(cdfCoef + y);
        return kSq * y;
    }
    y = std::exp(cdfCoef - k * y);
    return 1.0 - y;
}

double GeometricStableRand::pdfByLevy(double x) const
{
    if (mu == 0) {
        /// Invert parameter for case of -Levy
        if (beta < 0) {
            x = -x;
        }
        double halfSigmaInv = 0.5 / sigma;
        double z = halfSigmaInv * x;
        double sqrtZ = std::sqrt(z);
        double y = -std::erfc(sqrtZ);
        y *= std::exp(z);
        y += M_1_SQRTPI / sqrtZ;
        y *= halfSigmaInv;
        return y;
    }

    /// If mu != 0
    double a = x / sigma, c = RandMath::sign(beta) * mu / sigma;
    double h = std::sqrt(1 - 2 * c);
    double K = 1 - c - h;
    K *= a / (c * c);
    K = std::exp(K);
    K /= 2 * sigma * c * h;

    if (mu < 0) {
        double p =  std::exp(a * h / (c * c));
        if (a < 0)
            return -2 * K * (1 + h) * p * p;
        double sqrtHalfA = std::sqrt(0.5 * a);
        double g = (1 + h) / c;
        g *= sqrtHalfA;
        g = std::erfc(-g);
        double q = (1 - h) / c;
        q *= sqrtHalfA;
        q = std::erfc(q);

        /// g can be too small while p is too big
        double y = -(g * p) * p;
        y *= (1 + h);
        y -= (1 - h) * q;
        y *= K;
        return y;
    }

    // TODO: for mu > 0 and beta == -1
    return 0.0;
}

double GeometricStableRand::f(double x) const
{
    /// Laplace case
    if (alpha == 2)
        return pdfLaplace(x);

    /// Cut zeros for half-infinite distribution
    if (alpha < 1 && std::fabs(beta) == 1) {
        if (x < 0 && mu >= 0 && beta > 0)
            return 0.0;
        if (x > 0 && mu <= 0 && beta < 0)
            return 0.0;
    }

    /// Levy case
    if (alpha == 0.5 && std::fabs(beta) == 1)
        return pdfByLevy(x);

    if (x == 0) {
        if (alpha == 1)
            return INFINITY;
        if (mu == 0) {
            double coef = std::tgamma(1.0 - alphaInv);
            if (!std::isfinite(coef))
                return INFINITY;
            return Z.f(0) * coef / sigma;
        }

        // can be hashed
        double c = std::lgamma(alpha);
        double signMu = RandMath::sign(mu);
        c += (alpha - 1) * std::log(sigma);
        c -= (alpha + 1) * std::log(signMu * mu / sigma);
        c = std::exp(c);
        c *= std::sin(M_PI_2 * alpha);
        c /= M_PI;
        c *= alpha * (1 - signMu * beta);

        double int1 = RandMath::integral([this, c] (double z)
        {
            if (z <= 0 || z >= 1)
                return 0.0;
            double denominator = 1.0 - z;
            double t = z / denominator;
            double temp = sigma * std::pow(t, alphaInv);
            double par = -mu * t / temp;
            double y = Z.f(par) / temp;
            y -= c / std::pow(t, alpha);
            y *= std::exp(-t);
            return y / (denominator * denominator);
        },
        0, 1);

        double int2 = c * std::tgamma(1 - alpha); // hash
        return int1 + int2;
    }

    return RandMath::integral([this, x] (double z)
    {
        if (z <= 0 || z >= 1)
            return 0.0;
        double denominator = 1.0 - z;
        double t = z / denominator;
        double temp = sigma * std::pow(t, alphaInv);
        double par = (x - mu * t) / temp;
        double y = std::exp(-t) / temp * Z.f(par);
        return y / (denominator * denominator);
    },
    0, 1);
}

double GeometricStableRand::F(double x) const
{
    if (alpha == 2)
        return cdfLaplace(x);
    return RandMath::integral([this, x] (double z)
    {
        if (z <= 0)
            return 1.0;
        if (z >= 1) {
            if (mu == 0)
                return Z.F(0);
            if (alpha > 1)
                return (mu > 0) ? 0.0 : 1.0;
            double par = (alpha < 1) ? 0 : -mu / sigma;
            return Z.F(par);
        }
        double denominator = 1.0 - z;
        double t = z / denominator;
        double par = x - mu * t;
        par /= std::pow(t, alphaInv) * sigma;
        double y = std::exp(-t) * Z.F(par);
        return y / (denominator * denominator);
    },
    0, 1);
}

double GeometricStableRand::variateForUnityExponent() const
{
    double U = UniformRand::Variate(-M_PI_2, M_PI_2);
    double W1 = ExponentialRand::StandardVariate();
    double W2 = ExponentialRand::StandardVariate();
    double pi_2BetaU = M_PI_2 + beta * U;
    double X = logSigma;
    X += std::log(sigma * W2 * pi_2BetaU / (M_PI_2 * W1 * std::cos(U)));
    X *= beta;
    X += pi_2BetaU * std::tan(U);
    X *= M_2_PI * sigma;
    X += mu;
    X *= W2;
    return X;
}

double GeometricStableRand::variateForCommonExponent() const
{
    double U = UniformRand::Variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::StandardVariate();
    double alphaUB = alpha * U + B;
    double X = std::sin(alphaUB);
    double W_adj = W / std::cos(U - alphaUB);
    X *= W_adj;
    double Y = ExponentialRand::StandardVariate();
    W_adj *= std::cos(U);
    W_adj = Y / W_adj;
    double R = S + alphaInv * std::log(W_adj);
    X *= std::exp(R);
    return X * sigma + mu * Y;
}

double GeometricStableRand::variateByLevy(bool positive) const
{
    double W = ExponentialRand::StandardVariate();
    double X = LevyRand::StandardVariate();
    if (!positive)
        X = -X;
    X *= sigma * W;
    X += mu;
    return X * W;
}

double GeometricStableRand::variateByCauchy() const
{
    double W = ExponentialRand::StandardVariate();
    double X = CauchyRand::Variate(mu, sigma);
    return X * W;
}

double GeometricStableRand::Variate() const
{
    if (alpha == 2)
        return (mu == 0) ? LaplaceRand::Variate(0, sigma) : LaplaceRand::Variate(0, sigma, k);
    if (alpha == 0.5) {
        if (beta == 1)
            return variateByLevy(true);
        if (beta == -1)
            return variateByLevy(false);
    }
    if (alpha == 1)
        return (beta == 0) ? variateByCauchy() : variateForUnityExponent();
    return variateForCommonExponent();
}

void GeometricStableRand::Sample(std::vector<double> &outputData) const
{
    if (alpha == 2) {
        if (mu == 0) {
            for (double &var : outputData)
                var = LaplaceRand::Variate(0, sigma);
        }
        else {
            for (double &var : outputData)
                var = LaplaceRand::Variate(0, sigma, k);
        }
    }
    else if (alpha == 0.5 && std::fabs(beta) == 1) {
        if (beta > 0) {
            for (double & var : outputData)
                var = variateByLevy(true);
        }
        else {
            for (double & var : outputData)
                var = variateByLevy(false);
        }
    }
    else if (alpha == 1) {
        if (beta == 0) {
            for (double &var : outputData)
                var = variateByCauchy();
        }
        else {
            for (double &var : outputData)
                var = variateForUnityExponent();
        }
    }
    else {
        for (double &var : outputData)
            var = variateForCommonExponent();
    }
}

double GeometricStableRand::Variance() const
{
    if (alpha == 2) {
        double y = sigma * sigma / kSq;
        return (1.0 + kSq * kSq) * y;
    }
    return INFINITY;
}

std::complex<double> GeometricStableRand::CF(double t) const
{
    return (t == 0) ? 1.0 : 1.0 / (1.0 + psi(t));
}

double GeometricStableRand::Median() const
{
    if (alpha == 2) {
        double y = 0.5 * (1.0 / kSq + 1.0);
        y = std::log(y);
        return sigma * k * y;
    }
    return ContinuousDistribution::Median();
}

double GeometricStableRand::Mode() const
{
    if (alpha == 1 || alpha == 2)
        return 0.0;
    return ContinuousDistribution::Mode();
}

double GeometricStableRand::Skewness() const
{
    if (alpha == 2) {
        double kSq = k * k;
        double k4 = kSq * kSq;
        double k6 = k4 * kSq;
        double z = (k4 + 1);
        z *= std::sqrt(z);
        double y = 2 * (1 - k6);
        return y / z;
    }
    return NAN;
}

double GeometricStableRand::ExcessKurtosis() const
{
    if (alpha == 2)
    {
        if (k == 1)
            return 3.0;
        return NAN; // TODO!
    }
    return NAN;
}
