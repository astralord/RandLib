#include "RandomVariable.h"
#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision

RandomVariable::RandomVariable()
{
}

std::string RandomVariable::toStringWithPrecision(const double a_value, const int n)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

void RandomVariable::sample(QVector<double> &outputData)
{
    for (double &var : outputData)
        var = variate();
}

void RandomVariable::cdf(const QVector<double> &x, QVector<double> &y)
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = F(x[i]);
}
