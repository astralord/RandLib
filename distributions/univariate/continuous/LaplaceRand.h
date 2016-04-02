#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ExponentialRand.h"

/**
 * @brief The LaplaceRand class
 */
class RANDLIBSHARED_EXPORT LaplaceRand : public ContinuousDistribution
{
    double mu, sigma, k;
    double sigmaInv; /// 1 / sigma
    double kInv, kSq; /// 1 / k and k * k
    double pdfCoef; /// 1 / (sigma * (k + 1 / k))
    double cdfCoef; /// 1 / (1 + k * k)

public:
    LaplaceRand(double location = 0, double scale = 1, double asymmetry = 1);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    void setAsymmetry(double asymmetry);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return sigma; }
    inline double getAsymmetry() const { return k; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double scale);
    static double variate(double location, double scale, double asymmetry);

    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;

    /// Maximum likelihood estimation
    bool fitLocationMLE(const std::vector<double> &sample);
    bool fitScaleMLE(const std::vector<double> &sample);
    bool fitLocationAndScaleMLE(const std::vector<double> &sample);
    
    /// Method of moments
    bool fitLocationMM(const std::vector<double> &sample);
    bool fitScaleMM(const std::vector<double> &sample);
    bool fitLocationAndScaleMM(const std::vector<double> &sample);
};

#endif // LAPLACERAND_H
