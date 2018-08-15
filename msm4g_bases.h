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



#endif /* msm4g_bases_h */
