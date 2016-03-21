#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GammaRand class
 * Gamma distribution
 * X ~ \Gamma(k, \theta)
 */
class RANDLIBSHARED_EXPORT GammaRand : public ContinuousDistribution
{
protected:
    double alpha, theta, beta;
    double alphaInv; /// 1.0 / alpha
    double variateCoef; /// (e + alpha) / (alpha * e)
    double cdfCoef; /// 1.0 / gamma(alpha)
    double pdfCoef; /// 1.0 / (gamma(alpha) * theta ^ alpha)

private:
    double m, s, s_2, d, b, w, v, c; /// constants for sampling
    void setConstantsForGenerator();

public:
    GammaRand(double shape = 1, double scale = 1);
    virtual ~GammaRand() {}
    std::string name() override;

    void setParameters(double shape, double scale);
    inline double getShape() const { return alpha; }
    inline double getScale() const { return theta; }
    inline double getRate() const { return beta; }

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

    /// Maximum-likelihood estimation
    bool fitScaleMLE(const QVector<double> &sample);
    bool fitShapeAndScaleMLE(const QVector<double> &sample);
    
    /// Method of moments
    bool fitShapeMM(const QVector<double> &sample);
    bool fitScaleMM(const QVector<double> &sample);
    bool fitShapeAndScaleMM(const QVector<double> &sample);

    bool fitRateBayes(const QVector<double> &sample, GammaRand &priorDistribution);

    /**
     * @brief getInverseGammaFunction
     * @return 1.0 / Gamma(k)
     */
    inline double getInverseGammaFunction() { return cdfCoef; }
};


/**
 * @brief The ChiSquaredRand class
 * Chi-squared distribution
 * X ~ \chi^2(n)
 *
 * X ~ Gamma(0.5 * n, 2)
 */
class RANDLIBSHARED_EXPORT ChiSquaredRand : public GammaRand
{
public:
    explicit ChiSquaredRand(int degree = 1);
    std::string name() override;

    void setDegree(int degree);
    inline int getDegree() const;

protected:
    /// prohibit to use gamma's public getters and setters
    using GammaRand::setParameters;
    using GammaRand::getShape;
    using GammaRand::getScale;
};


/**
 * @brief The ErlangRand class
 * Erlang distibution
 * X ~ Erlang(k, l)
 *
 * X ~ Y_1 + Y_2 + ... + Y_k, where Y_i ~ Exp(l)
 * X ~ Gamma(k, 1/l)
 */
class RANDLIBSHARED_EXPORT ErlangRand : public GammaRand
{
public:
    ErlangRand(int shape = 1, double rate = 1);
    std::string name() override;

    void setParameters(int shape, double rate);
    inline int getShape() const;
    inline double getRate() const;

protected:
    /// prohibit to use gamma's public getters and setters
    using GammaRand::setParameters;
    using GammaRand::getShape;
    using GammaRand::getScale;
};


#endif // GAMMARAND_H
