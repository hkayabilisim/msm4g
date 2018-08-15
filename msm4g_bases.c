#include "msm4g_bases.h"

static double cubic_hermite0(double t) { return 0.5 * (1 - t) * (2 - t) * (2 - t); }
static double cubic_hermite1(double t) { return  (1 - t) * (1 + t - 1.5f * t * t); }
static double cubic_hermite2(double t) { return  (1 + t) * (1 - t - 1.5f * t * t); }
static double cubic_hermite3(double t) { return 0.5 * (1 + t) * (2 + t) * (2 + t); }

static double cubic_bspline0(double t) { return  (1./6) * (2-t) * (2-t) * (2-t);}
static double cubic_bspline1(double t) { return  (2./3) + t*t*(-1 + 0.5*t);     }
static double cubic_bspline2(double t) { return  (2./3) + t*t*(-1 - 0.5*t);     }
static double cubic_bspline3(double t) { return  (1./6) * (2+t) * (2+t) * (2+t);}

static const struct BaseFunction _CubicHermite = { 4, cubic_hermite0, cubic_hermite1, cubic_hermite2, cubic_hermite3 };
static const struct BaseFunction _CubicBSpline = { 4, cubic_bspline0, cubic_bspline1, cubic_bspline2, cubic_bspline3 };

const BaseFunction *CubicHermite = & _CubicHermite;
const BaseFunction *CubicBSpline = & _CubicBSpline;


