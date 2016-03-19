#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GammaRand class
 */
class RANDLIBSHARED_EXPORT GammaRand : public ContinuousDistribution
{
protected:
    double k, theta, thetaInv;
    
private:
    double kInv; /// 1.0 / k
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

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool checkValidity(const QVector<double> &sample);

    bool fitScale_MLE(const QVector<double> &sample);
    virtual bool fit_MLE(const QVector<double> &sample);
    
    bool fitShape_MM(const QVector<double> &sample);
    bool fitScale_MM(const QVector<double> &sample);
    virtual bool fit_MM(const QVector<double> &sample);

    /**
     * @brief getInverseGammaFunction
     * @return 1.0 / Gamma(k)
     */
    inline double getInverseGammaFunction() { return cdfCoef; }
};

#endif // GAMMARAND_H
