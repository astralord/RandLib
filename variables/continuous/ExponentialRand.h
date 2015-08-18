#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The ExponentialRand class
 */
class RANDLIBSHARED_EXPORT ExponentialRand : public ContinuousRand
{
    double lambda, beta;

    //TODO: find a way to initialize them without dummy
    /// Tables for ziggurat
    static double stairWidth[257], stairHeight[256];
    static double constexpr x1 = 7.69711747013104972;
    static bool dummy;
    static bool setupTables();

public:
    ExponentialRand(double rate = 1);
    virtual void setName() override;

    void setRate(double rate);
    inline double getRate() const { return lambda; }
    inline double getScale() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    std::complex<double> CF(double t) const override;
    double variate() const override;

    static double variate(double lambda);
    static double standardVariate();

    double E() const override { return beta; }
    double Var() const override { return beta * beta; }

    static constexpr double Skewness() { return 2.0; }
    static constexpr double ExcessiveKurtosis() { return 6.0; }

    inline double Entropy() const { return 1 - std::log(lambda); }

    bool fitToData(const QVector<double> &sample);
};

#endif // EXPONENTIALRAND_H
