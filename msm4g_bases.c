#include "msm4g_bases.h"

double msm4g_bases_bspline(int k, double u) {
    if (k == 1) {
        if (u >= 0 && u < 1)
            return 1.0;
        else
            return 0.0;
    } else {
        return u / (double) (k - 1) * msm4g_bases_bspline(k - 1, u)
        + (k - u) / (double) (k - 1) * msm4g_bases_bspline(k - 1, u - 1.0);
    }
}

/* From Section III.B.1 of Bob's B-spline article */
double msm4g_bases_bsplineprime(int k, double u) {
    return msm4g_bases_bspline(k - 1, u) - msm4g_bases_bspline(k - 1, u - 1);
}
