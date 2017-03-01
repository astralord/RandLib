#ifndef GAMMARAND_H
#define GAMMARAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GammaRand class
 * Gamma distribution
 *
 * Notation X ~ Gamma(α, β)
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
    double logBeta; /// log(β)

private:
    double t, b, alphaInv; /// constants for sampling
    void SetConstantsForGenerator();

public:
    GammaRand(double shape = 1, double rate = 1);
    virtual ~GammaRand() {}

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }
    void SetParameters(double shape, double rate);
    inline double GetShape() const { return alpha; }
    inline double GetScale() const { return theta; }
    inline double GetRate() const { return beta; }

    double f(double x) const override;
    double logf(double x) const override;
    double F(double x) const override;
    double S(double x) const override;
    
private:

    enum GENERATOR_ID {
        INTEGER_SHAPE, /// Erlang distribution for α = 1, 2, 3
        ONE_AND_A_HALF_SHAPE, /// α = 1.5
        SMALL_SHAPE, /// α < 0.34
        FISHMAN, /// 1 < α < 1.2
        MARSAGLIA_TSANG /// 0.34 < α < 1 or α >= 1.2
    };

    static GENERATOR_ID GetIdOfUsedGenerator(double shape);

    static double variateThroughExponentialSum(int shape);
    static double variateForShapeOneAndAHalf();
    double variateBest() const;
    static double variateAhrensDieter(double shape);
    static double variateFishman(double shape);
    static double variateMarsagliaTsang(double shape);
    
public:
    static double StandardVariate(double shape);
    static double Variate(double shape, double rate);

    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double initRootForSmallP(double r) const;
    double initRootForLargeP(double logQ) const;
    double initRootForSmallShape(double p) const;
    double initRootForLargeShape(double p) const;
    double df(double x) const;
    double quantileInitialGuess(double p) const;
    double quantileInitialGuess1m(double p) const;
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    /// Maximum-likelihood estimation
    void FitScaleMLE(const std::vector<double> &sample);
    void FitMLE(const std::vector<double> &sample);
    
    /// Method of moments
    void FitShapeMM(const std::vector<double> &sample);
    void FitScaleMM(const std::vector<double> &sample);
    void FitMM(const std::vector<double> &sample);

    /// Bayes estimation
    GammaRand FitRateBayes(const std::vector<double> &sample, const GammaRand &priorDistribution);

    /**
     * @brief GetLogGammaFunction
     * @return log(Gamma(α))
     */
    inline double GetLogGammaFunction() const { return -mLgammaShape; }
    /**
     * @brief GetLogRate
     * @return log(β)
     */
    inline double GetLogRate() const { return logBeta; }
};


/**
 * @brief The ChiSquaredRand class
 * Chi-squared distribution
 *
 * Notation: X ~ Chi^2(n)
 *
 * Related distributions:
 * X ~ Gamma(0.5 * n, 0.5)
 */
class RANDLIBSHARED_EXPORT ChiSquaredRand : public GammaRand
{
public:
    explicit ChiSquaredRand(int degree = 1);
    std::string Name() const override;

    void SetDegree(int degree);
    inline int GetDegree() const { return static_cast<int>(alpha + alpha); }

protected:
    /// prohibit to use gamma's public Getters and Setters
    using GammaRand::SetParameters;
    using GammaRand::GetShape;
    using GammaRand::GetScale;
};


/**
 * @brief The ErlangRand class
 * Erlang distibution
 *
 * Notation: X ~ Erlang(k, l)
 *
 * Related distributions:
 * X ~ Y_1 + Y_2 + ... + Y_k, where Y_i ~ Exp(l)
 * X ~ Gamma(k, l)
 */
class RANDLIBSHARED_EXPORT ErlangRand : public GammaRand
{
public:
    ErlangRand(int shape = 1, double rate = 1);
    std::string Name() const override;

    inline int GetShape() const;
    inline double GetRate() const;

protected:
    /// prohibit to use gamma's public Getters and Setters
    using GammaRand::SetParameters;
    using GammaRand::GetShape;
    using GammaRand::GetScale;
};


#endif // GAMMARAND_H
