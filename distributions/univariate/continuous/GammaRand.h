#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GammaRand class
 * Gamma distribution
 *
 * Notation X ~ Gamma(k, β)
 *
 * Related distributions:
 * If X ~ Gamma(1, β), then X ~ Exp(β)
 */
class RANDLIBSHARED_EXPORT GammaRand : public ContinuousDistribution
{
protected:
    double alpha, theta, beta;
    double mLgammaShape; /// -lgamma(α)
    double pdfCoef; ///  -lgamma(α) - α * log(θ)

private:
    double t, b, alphaInv; /// constants for sampling
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

    enum GENERATOR_ID {
        INTEGER_SHAPE, /// Erlang distribution for α = 1, 2, 3
        ONE_AND_A_HALF_SHAPE, /// α = 1.5
        SMALL_SHAPE, /// α < 0.34
        FISHMAN, /// 1 < α < 1.2
        MARSAGLIA_TSANG /// 0.34 < α < 1 or α >= 1.2
    };

    static GENERATOR_ID getIdOfUsedGenerator(double shape);

    static double variateThroughExponentialSum(int shape);
    static double variateForShapeOneAndAHalf();
    double variateBest() const;
    static double variateAhrensDieter(double shape);
    static double variateFishman(double shape);
    static double variateMarsagliaTsang(double shape);
    
public:
    static double standardVariate(double shape);
    static double variate(double shape, double rate);

    double variate() const override; 
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double QuantileImpl(double p) const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    /// Maximum-likelihood estimation
    bool fitScaleMLE(const std::vector<double> &sample);
    bool fitMLE(const std::vector<double> &sample);
    
    /// Method of moments
    bool fitShapeMM(const std::vector<double> &sample);
    bool fitScaleMM(const std::vector<double> &sample);
    bool fitMM(const std::vector<double> &sample);

    /// Bayes estimation
    bool fitRateBayes(const std::vector<double> &sample, GammaRand &priorDistribution);

    /**
     * @brief getLogGammaFunction
     * @return log(Gamma(α))
     */
    inline double getLogGammaFunction() const { return -mLgammaShape; }
};


/**
 * @brief The ChiSquaredRand class
 * Chi-squared distribution
 *
 * Notation: X ~ Chi^2(n)
 *
 * Related distributions:
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
 *
 * Notation: X ~ Erlang(k, l)
 *
 * Related distributions:
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
