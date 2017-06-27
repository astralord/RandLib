#ifndef NUMERICMATH
#define NUMERICMATH

#include "RandMath.h"

/// Numerical procedures

namespace RandMath
{

/**
 * @fn integral
 * @param funPtr integrand
 * @param a lower boundary
 * @param b upper boundary
 * @param epsilon tolerance
 * @param maxRecursionDepth how deep should the algorithm go
 * @return
 */
double integral(const std::function<double (double)> &funPtr, double a, double b,
                            double epsilon = 1e-11, int maxRecursionDepth = 11);

/**
 * @fn findRoot
 * Newton's root-finding procedure,
 * using first and second derivatives
 * @param funPtr mapping x |-> (f(x), f'(x), f''(x))
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoubleTriplet (double)> &funPtr, double & root, double funTol = 1e-10, double stepTol = 1e-6);

/**
 * @fn findRoot
 * Newton's root-finding procedure,
 * using first derivative
 * @param funPtr mapping x |-> (f(x), f'(x))
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoublePair (double)> &funPtr, double & root, double funTol = 1e-10, double stepTol = 1e-6);

/**
 * @fn findRoot
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
 * @fn findMin
 * Combined Brent's method
 * @param funPtr
 * @param abc lower boundary / middle / upper boundary
 * @param fabc funPtr(abc)
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return true if success
 */
bool findMin(const std::function<double (double)> &funPtr, const DoubleTriplet &abc, const DoubleTriplet &fabc, double &root, double epsilon = 1e-8);

/**
 * @fn findMin
 * Combined Brent's method
 * @param funPtr
 * @param closePoint point that is nearby minimum
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return true if success
 */
bool findMin(const std::function<double (double)> &funPtr, double closePoint, double &root, double epsilon = 1e-8);

}

#endif // NUMERICMATH

