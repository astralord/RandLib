#include "NormalRand.h"

unsigned long NormalRand::kn[128] = {0};
double NormalRand::wn[128] = {0};
double NormalRand::fn[128] = {0};
const bool NormalRand::dummy = NormalRand::setupTables();

NormalRand::NormalRand(double mean, double var)
{
    setMean(mean);
    setVar(var);
}

bool NormalRand::setupTables()
{
    const double m1 = 2147483648.0;
    double dn = 3.442619855899, tn = dn, vn = 9.91256303526217e-3;

    fn[0] = 1.;
    fn[127] = std::exp(-.5 * dn * dn);

    double q = vn / fn[127];

    kn[0] = (dn / q) * m1;
    kn[1] = 0;

    wn[0] = q / m1;
    wn[127] = dn / m1;

    for (size_t i = 126; i >= 1; --i)
    {
        dn = -std::log(vn / dn + std::exp(-.5 * dn * dn));
        dn = std::sqrt(dn + dn);
        kn[i + 1] = (dn / tn) * m1;
        tn = dn;
        fn[i] = std::exp(-.5 * dn * dn);
        wn[i] = dn / m1;
    }
    return true;
}

void NormalRand::setMean(double mean)
{
    mu = mean;
}

void NormalRand::setSigma(double rootVar)
{
    sigma = std::max(rootVar, MIN_POSITIVE);
    sigmaSqrt2Inv = M_SQRT1_2 / sigma;
}

double NormalRand::f(double x) const
{
    double y = x - mu; /// x - mu
    y *= sigmaSqrt2Inv; /// (x - mu) / (sigma * sqrt(2))
    y *= y; /// (((x - mu) / sigma) ^ 2) / 2
    y = std::exp(-y); /// exp((((x - mu) / sigma) ^ 2) / 2)
    return M_1_SQRTPI * sigmaSqrt2Inv * y; /// exp((((x - mu) / sigma) ^ 2) / 2) / (sigma * sqrt(2pi))
}

double NormalRand::F(double x) const
{
    double y = x - mu; /// x - mu
    y *= sigmaSqrt2Inv; /// (x - mu) / (sigma * sqrt(2))
    y = std::erf(y); /// erf((x - mu) / (sigma * sqrt(2)))
    ++y; /// 1 + erf((x - mu) / (sigma * sqrt(2)))
    return .5 * y; /// (1 + erf((x - mu) / (sigma * sqrt(2)))) / 2
}

double NormalRand::ziggurat()
{
    int iter = 0;
    do {
        long hz = BasicRandGenerator::getRand();
        unsigned long iz = hz & 127;
        double x = hz * wn[iz];
        if (std::fabs(hz) < kn[iz])
            return x;

        if (iz == 0) /// Handle the base strip
        {
            double y;
            do {
                y = -std::log(U.variate());
                x = -std::log(U.variate()) * 0.2904764; /// 1.0 / 3.44262
            } while (y + y < x * x);

            //TODO: maybe we need more accuracy?
            x += 3.44262; /// + start of the right tail
            return (hz > 0) ? x : -x;
        }

        /// Handle the wedges of other strips
        if (fn[iz] + U.variate() * (fn[iz - 1] - fn[iz]) < std::exp(-.5 * x * x))
            return x;
    } while (++iter <= 1e9); /// one billion should be enough
    return 0;
}

double NormalRand::variate()
{
    return mu + sigma * ziggurat();
}

bool NormalRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate mu
    double average = 0.0;
    for (double var : sample) {
        average += var;
    }
    average /= sample.size();

    /// Calculate sigma
    double deviation = 0.0;
    for (double var : sample) {
        double currDev = (var - average);
        deviation += currDev * currDev;
    }
    deviation /= std::max(sample.size() - 1, 1);

    setMean(average);
    setSigma(std::sqrt(deviation));
    return true;
}
