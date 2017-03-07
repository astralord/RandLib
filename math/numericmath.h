#ifndef NUMERICMATH
#define NUMERICMATH

#include "RandMath.h"

/// Numerical procedures

namespace RandMath
{

/**
 * @brief integral
 * @param funPtr integrand
 * @param a lower boundary
 * @param b upper boundary
 * @param epsilon tolerance
 * @param maxRecursionDepth how deep should the algorithm go
 * @return
 */
long double integral(const std::function<double (double)> &funPtr, double a, double b,
                            double epsilon = 1e-11, int maxRecursionDepth = 11);

/**
 * @brief findRoot
 * Newton's root-finding procedure,
 * using first and second derivatives
 * @param funPtr mapping x |-> (f(x), f'(x), f''(x))
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoubleTriplet (double)> &funPtr, double & root, double funTol = 1e-10, double stepTol = 1e-6);

/**
 * @brief findRoot
 * Newton's root-finding procedure,
 * using first derivative
 * @param funPtr mapping x |-> (f(x), f'(x))
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoublePair (double)> &funPtr, double & root, double funTol = 1e-10, double stepTol = 1e-6);

/**
 * @brief findRoot
 * Brent's root-finding procedure
 * @param funPtr mapping x |-> f(x)
 * @param a lower boundary
 * @param b upper boundary
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<double (double)> &funPtr, double a, double b, double & root, double epsilon = 1e-8);

/**
 * @brief findMin
 * Combined Brent's method
 * @param funPtr
 * @param a lower boundary
 * @param c upper boundary
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return
 */
bool findMin(const std::function<double (double)> &funPtr, const DoubleTriplet &abc, const DoubleTriplet &fabc, double &root, double epsilon = 1e-8);

/**
 * @brief findMin
 * Combined Brent's method
 * @param funPtr
 * @param closePoint point that is nearby minimum
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return
 */
bool findMin(const std::function<double (double)> &funPtr, double closePoint, double &root, double epsilon = 1e-8);

}

#endif // NUMERICMATH

