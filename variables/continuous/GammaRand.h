#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

#include <functional>

/**
 * @brief The GammaRand class
 */
class RANDLIBSHARED_EXPORT GammaRand : public ContinuousRand
{
protected:
    double k, theta;
    
private:
    double kInv, thetaInv; /// 1.0 / k and 1.0 / theta
    double variateCoef; /// (e + k) / (k * e)
    double cdfCoef; /// 1.0 / gamma(k)
    double pdfCoef; /// 1.0 / (gamma(k) * theta ^ k)
    double m, s, s_2, d, b, w, v, c; /// constants for sampling
    void setConstantsForGenerator();

public:
    GammaRand(double shape = 1, double scale = 1);
    virtual ~GammaRand() {}
    std::string name() override;

    void setParameters(double shape, double scale);
    inline double getShape() const { return k; }
    inline double getScale() const { return theta; }

    double f(double x) const override;
    double F(double x) const override;
    
private:
    double variateForIntegerShape() const;     /// Erlang distribution (use for k < 5)
    double variateForHalfIntegerShape() const; /// GA algorithm for k = [n] + 0.5 and k < 5
    double variateForSmallShape() const;       /// GS algorithm for small k < 1
    double variateForMediumShape() const;      /// GP algorithm for 1 < k < 3
    double variateForLargeShape() const;       /// GO algorithm for the most common case k > 3
    
public:
    double variate() const override;

    void sample(QVector<double> &outputData) const override;

    double Mean() const override { return k * theta; }
    double Variance() const override { return k * theta * theta; }

    std::complex<double> CF(double t) const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    /**
     * @brief getInverseGammaFunction
     * @return 1.0 / Gamma(k)
     */
    inline double getInverseGammaFunction() { return cdfCoef; }
};

#endif // GAMMARAND_H
