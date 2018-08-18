#ifndef msm4g_bases_h
#define msm4g_bases_h

#include <stdio.h>

#define MAX_POLY_DEGREE 20

typedef struct BaseFunction
{
    const int p;
    double (*region[MAX_POLY_DEGREE]) (double x);
} BaseFunction;

extern const BaseFunction *CubicHermite;
extern const BaseFunction *CubicBSpline;

/** @brief Evaluates B-Spline bases
 *
 * Calculates the B-Spline bases of order k at any given 
 * point u. It uses the recursive definition given in
 * Eq 7 in "The Journal of Chemical Physics 144, 114112 (2016); doi: 10.1063/1.4943868".
 * 
 * @todo Using the recursive definition is easy and appearently not too costly
 * but it is one of the first places to look for if overall MSM has some performance
 * issues.
 *
 * @param[in] k Order of B-Spline
 * @param[in] u Evaluation point
 * @return B-Spline of order k evaluated at point u
 */
double msm4g_bases_bspline(int k,double u);

/** @brief Evaluates the derivative B-Spline bases
 *
 * Calculates the derivative of B-Spline bases of order k at any given
 * point u. It uses the recursive definition given in
 * Section II.B.1 in "The Journal of Chemical Physics 144, 114112 (2016); doi: 10.1063/1.4943868".
 *
 * @todo Like its non-derivative version, it also uses a recursive definition.
 * So it is also one of the first places to look for if overall MSM has some performance
 * issues.
 *
 * @param[in] k Order of B-Spline
 * @param[in] u Evaluation point
 * @return Derivative of B-Spline of order k evaluated at point u
 */
double msm4g_bases_bsplineprime(int k,double u);

#endif /* msm4g_bases_h */
