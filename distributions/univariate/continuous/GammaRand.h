#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GammaRand class
 * Gamma distribution
 * X ~ \Gamma(k, θ)
 */
class RANDLIBSHARED_EXPORT GammaRand : public ContinuousDistribution
{
protected:
    double alpha, theta, beta;
    double cdfCoef; /// -lgamma(α)
    double pdfCoef; ///  -lgamma(α) - α * log(θ)

private:
    double m, s, s_2, d, b, w, v, c; /// constants for sampling
    void setConstantsForGenerator();

public:
    GammaRand(double shape = 1, double rate = 1);
    virtual ~GammaRand() {}

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }
    void setParameters(double shape, double rate);
    inline double getShape() const { return alpha; }
    inline double getScale() const { return theta; }
    inline double getRate() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    
private:
    static double variateForIntegerShape(int shape);     /// Erlang distribution (use for k < 5)
    static double variateForHalfIntegerShape(int shape); /// GA algorithm for k = [n] + 0.5 and k < 5
    static double variateForSmallShape(double shape);    /// GS algorithm for small k < 1
    static double variateForMediumShape(double shape);   /// GP algorithm for 1 < k < 3
    double variateForLargeShape() const;                 /// GO algorithm for the most common case k > 3
    static double variateForLargeShape(double shape);
    
public:
    static double standardVariate(double shape);
    static double variate(double shape, double rate);

    double variate() const override; 
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool checkValidity(const std::vector<double> &sample);

    /// Maximum-likelihood estimation
    bool fitScaleMLE(const std::vector<double> &sample);
    bool fitShapeAndScaleMLE(const std::vector<double> &sample);
    
    /// Method of moments
    bool fitShapeMM(const std::vector<double> &sample);
    bool fitScaleMM(const std::vector<double> &sample);
    bool fitShapeAndScaleMM(const std::vector<double> &sample);

    /// Bayes estimation
    bool fitRateBayes(const std::vector<double> &sample, GammaRand &priorDistribution);

    /**
     * @brief getLogGammaFunction
     * @return 1.0 / Gamma(k)
     */
    inline double getLogGammaFunction() const { return -cdfCoef; }
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
    std::string name() const override;

    void setDegree(int degree);
    inline int getDegree() const { return static_cast<int>(alpha + alpha); }

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
    std::string name() const override;

    inline int getShape() const;
    inline double getRate() const;

protected:
    /// prohibit to use gamma's public getters and setters
    using GammaRand::setParameters;
    using GammaRand::getShape;
    using GammaRand::getScale;
};


#endif // GAMMARAND_H
