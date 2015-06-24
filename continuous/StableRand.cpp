#include "StableRand.h"
#include <QDebug>

StableRand::StableRand(double exponent, double skewness, double scale, double location) :    
    U(-M_PI_2, M_PI_2), N(0, M_SQRT2)
{
    setParameters(exponent, skewness, scale, location);
}

void StableRand::setParameters(double exponent, double skewness, double scale, double location)
{
    alpha = std::min(exponent, 2.0);
    alpha = std::max(alpha, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
    alpha_alpham1 = alpha / (alpha - 1.0);
    alpham1Inv = alpha_alpham1 - 1.0;

    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);

    sigma = std::max(scale, MIN_POSITIVE);
    mu = location;

    /// Should be cautious, known distributions in priority
    if (std::fabs(alpha - 2) < MIN_POSITIVE)
        alpha = 2;
    else if (std::fabs(alpha - 1) < MIN_POSITIVE)
        alpha = 1;
    else if (std::fabs(alpha - .5) < MIN_POSITIVE)
        alpha = .5;

    if (std::fabs(beta) < MIN_POSITIVE)
        beta = 0;
    else if (std::fabs(beta - 1) < MIN_POSITIVE)
        beta = 1;
    else if (std::fabs(beta + 1) < MIN_POSITIVE)
        beta = -1;

    if (alpha == 2) /// X ~ Normal(mu, 2sigma^2)
    {
        valuePtr = std::bind(&StableRand::valueNormal, this);
        pdfPtr = std::bind(&StableRand::pdfNormal, this, std::placeholders::_1);
        cdfPtr = std::bind(&StableRand::cdfNormal, this, std::placeholders::_1);
    }
    else if (alpha == 1)
    {
        if (beta == 0) /// X ~ Cauchy(mu, sigma)
        {
            valuePtr = std::bind(&StableRand::valueCauchy, this);
            pdfPtr = std::bind(&StableRand::pdfCauchy, this, std::placeholders::_1);
            cdfPtr = std::bind(&StableRand::cdfCauchy, this, std::placeholders::_1);
        }
        else /// just alpha == 1
        {
            logSigma = std::log(sigma);
            theta0 = 0;
            A = 1;

            valuePtr = std::bind(&StableRand::valueForAlphaEqualOne, this);
            pdfPtr = std::bind(&StableRand::pdfForAlphaEqualOne, this, std::placeholders::_1);
            cdfPtr = std::bind(&StableRand::cdfForAlphaEqualOne, this, std::placeholders::_1);
        }
    }
    else if (alpha == .5)
    {
        if (beta == 1) /// X ~ Levy(mu, sigma)
        {
            valuePtr = std::bind(&StableRand::valueLevy, this);
            pdfPtr = std::bind(&StableRand::pdfLevy, this, std::placeholders::_1);
            cdfPtr = std::bind(&StableRand::cdfLevy, this, std::placeholders::_1);
        }
        else if (beta == -1) /// -X ~ Levy(mu, sigma)
        {
            valuePtr = std::bind(&StableRand::valueLevyNegative, this);
            pdfPtr = std::bind(&StableRand::pdfLevyNegative, this, std::placeholders::_1);
            cdfPtr = std::bind(&StableRand::cdfLevyNegative, this, std::placeholders::_1);
        }
    }
    else if (alpha != 1) /// Common case
    {
        B = beta * std::tan(M_PI_2 * alpha);
        zeta = -B;
        S = std::pow(1 + B * B, .5 * alphaInv);
        B = std::atan(B);
        pdfCoef = M_1_PI * alpha / (std::fabs(1 - alpha) * sigma);
        theta0 = alphaInv * B;
        A = std::pow(qFastCos(B), alpham1Inv);

        valuePtr = std::bind(&StableRand::valueForCommonAlpha, this);
        pdfPtr = std::bind(&StableRand::pdfForCommonAlpha, this, std::placeholders::_1);
        cdfPtr = std::bind(&StableRand::cdfForCommonAlpha, this, std::placeholders::_1);
    }
}

double StableRand::value()
{
    /// Get standard value
    double rv = 0;
    int iter = 0;
    do {
        rv = valuePtr();
        ++iter;
    } while ((std::isnan(rv) || std::isinf(rv)) &&
             iter < 10); /// if we got nan 10 times - we have a problem, get out
    /// and return shifted-scaled
    return mu + sigma * rv;
}

double StableRand::valueForCommonAlpha()
{
    double V = U.value();
    double W = E.value();
    double alphaVB = alpha * V + B;
    double rv = S * qFastSin(alphaVB); /// S * sin(alpha * V + B)
    double W_adj = W / qFastCos(V - alphaVB);
    rv *= W_adj; /// S * sin(alpha * V + B) * W / cos((1 - alpha) * V - B)
    rv *= std::pow(W_adj * qFastCos(V), -alphaInv);/// S * sin(alpha * V + B) * W / cos((1 - alpha) * V - B) /
                                                   /// ((W * cos(V) / cos((1 - alpha) * V - B)) ^ (1 / alpha))
    return rv;
}

double StableRand::valueForAlphaEqualOne()
{
    double V = U.value();
    double W = E.value();
    double pi_2BetaV = M_PI_2 + beta * V;
    double rv = pi_2BetaV * std::tan(V);
    rv -= std::log(W * qFastCos(V) / pi_2BetaV);
    rv += logSigma;
    rv *= M_2_PI * beta;
    return rv;
}

double StableRand::f(double x) const
{
    return pdfPtr(x);
}

double StableRand::pdfForCommonAlpha(double x)
{
    if (std::fabs(x - mu) < MIN_POSITIVE)
    {
        double numen = std::tgamma(1 + alphaInv);
        numen *= qFastCos(theta0);
        double denom = M_PI * std::pow(1 + theta0 * theta0, .5 * alphaInv);
        return numen / denom;
    }

    x = (x - mu) / sigma + zeta;
    double x_adj = std::pow(x - zeta, alpha_alpham1);

    if (x > zeta)
    {
        if (alpha < 1 && beta == -1)
            return 0;

        // NOT QUICK, BUT RIGHT
        /*g = @(theta) xshift(:) .* V(theta) - 1;
        R = repmat([-theta0, pi/2 ],numel(xshift),1);
        if abs(beta) < 1
            theta2 = bisectionSolver(g,R,alpha);
        else
            theta2 = bisectionSolver(g,R,alpha,beta,xshift);
        end
        theta2 = reshape(theta2,size(xshift));
        % change variables so the two integrals go from
        % 0 to 1/2 and 1/2 to 1.
        theta2shift1 = 2*(theta2 + theta0);
        theta2shift2 = 2*(pi/2 - theta2);
        g1 = @(theta)  xshift .* ...
            V(theta2shift1 * theta - theta0);
        g2 = @(theta)  xshift .* ...
            V(theta2shift2 * (theta - .5) + theta2);
        zexpz = @(z) max(0,z .* exp(-z)); % use max incase of NaN

        p(x > zeta) = c2 .* ...
            (theta2shift1 .* quadv(@(theta) zexpz( g1(theta) ),...
                                    0 , .5, tol) ...
           + theta2shift2 .* quadv(@(theta) zexpz( g2(theta) ),...
                                   .5 , 1, tol) );  */

        // QUICK
        double y = integralPDF(-theta0, MIN_POSITIVE, 10, x_adj);
        return pdfCoef / (x - zeta) * y;
    }
    else
    {
        if (alpha < 1 && beta == 1)
            return 0;
        double y = integralPDF(theta0, MIN_POSITIVE, 10, -x_adj);
        return -pdfCoef / (x - zeta) * y;
    }
}

double StableRand::integrand(double gamma, double gamma0, double coef) const
{
    double cosGammaInv = 1.0 / std::max(qFastCos(gamma), MIN_POSITIVE);
    double gammaAdj = alpha * (gamma - gamma0);
    double y = std::max(cosGammaInv * qFastSin(gammaAdj), MIN_POSITIVE);
    y = std::pow(y, -alpha_alpham1);
    y *= qFastCos(gammaAdj - gamma);
    double u = coef * A * y * cosGammaInv;
    return u * std::exp(-u);
}

double StableRand::adaptiveSimpsonsAux(double coef, double gamma0, double a, double b, double epsilon, double S,
                                       double fa, double fb, double fc, int bottom) const
{
    // TODO: rewrite recursion into loop
    double c = .5 * (a + b), h = (b - a) / 12.0;
    double d = .5 * (a + c), e = .5 * (c + b);
    double fd = integrand(d, gamma0, coef), fe = integrand(e, gamma0, coef);
    double Sleft = h * (fa + 4 * fd + fc);
    double Sright = h * (fc + 4 * fe + fb);

    double S2 = Sleft + Sright;
    double dev = (S2 - S) / 15.0;
    if (bottom <= 0 || std::fabs(dev) <= epsilon)
        return S2 + dev;

    epsilon *= .5;
    --bottom;

    return adaptiveSimpsonsAux(coef, gamma0, a, c, epsilon, Sleft, fa, fc, fd, bottom) +
           adaptiveSimpsonsAux(coef, gamma0, c, b, epsilon, Sright, fc, fb, fe, bottom);
}

double StableRand::integralPDF(double gamma0, double epsilon, int maxRecursionDepth, double coef) const
{
    double c = .5 * (gamma0 + M_PI_2), h = (M_PI_2 - gamma0) / 6.0;
    double fa = integrand(gamma0, gamma0, coef), fb = integrand(M_PI_2, gamma0, coef), fc = integrand(c, gamma0, coef);
    double S = h * (fa + 4 * fc + fb);
    return adaptiveSimpsonsAux(coef, gamma0, gamma0, M_PI_2, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

double StableRand::F(double x) const
{
    return cdfPtr(x);
}
