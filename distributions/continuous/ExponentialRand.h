#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The ExponentialRand class
 */
class RANDLIBSHARED_EXPORT ExponentialRand : public ContinuousDistribution
{
    double lambda, beta;

    //TODO: find a way to initialize them without dummy
    /// Tables for ziggurat
    static double stairWidth[257], stairHeight[256];
    static double constexpr x1 = 7.69711747013104972;
    static bool dummy;
    static bool setupTables();

public:
    explicit ExponentialRand(double rate = 1);
    std::string name() override;

    void setRate(double rate);
    inline double getRate() const { return lambda; }
    inline double getScale() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double lambda);
    static double standardVariate();

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;

    double Moment(int n) const;

    bool checkValidity(const QVector<double> &sample);

    bool fit_MM(const QVector<double> &sample);
    bool fit_MLE(const QVector<double> &sample);
    bool fit_UMVU(const QVector<double> &sample);
};

#endif // EXPONENTIALRAND_H
