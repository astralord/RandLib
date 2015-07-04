#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

#include <functional>

/**
 * @brief The GammaRand class
 */
class RANDLIBSHARED_EXPORT GammaRand : public ContinuousRand
{
    double k, theta;
    double kInv, thetaInv; /// 1.0 / k and 1.0 / theta
    double valueCoef; /// (e + k) / (k * e)
    double cdfCoef; /// 1.0 / gamma(k)
    double pdfCoef; /// 1.0 / (gamma(k) * theta ^ k)
    UniformRand U;
    ExponentialRand E;
    NormalRand N;

    double m, s, s_2, d, b, w, v, c; /// constants for sampling
    void setConstants();

public:
    GammaRand(double shape = 1, double scale = 1);

    void setParameters(double shape, double scale);
    inline double getShape() const { return k; }
    inline double getScale() const { return theta; }

    virtual double value() override;
    virtual double f(double x) const override;
    virtual double F(double x) const override;

    double M() const override { return k * theta; }
    double Var() const override { return k * theta * theta; }

private:
    std::function<double ()> valuePtr;
    double valueForIntegerShape();     /// Erlang distribution
    double valueForHalfIntegerShape(); /// GA algorithm for k = [n] + 0.5
    double valueForSmallShape();       /// GS algorithm for small k < 1
    double valueForMediumShape();      /// GP algorithm for 1 < k < 3
    double valueForLargeShape();       /// GO algorithm for the most common case k > 3
};

#endif // GAMMARAND_H
