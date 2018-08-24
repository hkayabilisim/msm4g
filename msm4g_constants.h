//
//  msm4g_constants.h
//  msm4g
//
//  Created by Huseyin Kaya on 3/29/17.
//  Copyright © 2017 Huseyin Kaya. All rights reserved.
//

#ifndef msm4g_constants_h
#define msm4g_constants_h


/** @brief 1D nesting stencil for Cubic B-Spline  */
//const double STENCIL_NESTING_1D_CUBIC[3]   = { 6./8 ,  4./8 , 1./8        };
/** @brief 1D nesting stencil for Quintic B-Spline  */
//const double STENCIL_NESTING_1D_QUINTIC[4] = {20./32, 15./32, 6./32, 1./32};

/** @brief Specifices the maximum number of omega values.
 * 
 * The relation between the row index (i)
 * of OMEGA array and p is: p = (i+2)*2 or i = p/2 - 2. So, the first row is for p=4, the second is for p=6, and so on.
 * You can access the omega values for p=4 (Cubic) and p=6 (Quintic) as follows:
 * @code
   for (int p=4 ; p <= 6; p = p + 2)
   {
     int i = p/2 - 2;
     for (int mu=0; mu < MAX_MU; mu++ )
        printf("OMEGA[%d][%d] (p=%d) : %25.16e\n",i,mu,p,OMEGA[i][mu]);
   }
 * @endcode
 */
#define MAX_MU           100
#define LMAX             20
#define MAX_POLY_DEGREE  20

#define FACE_MAXLEN 40000
#define TOL_DIRECT  1E-8
#define TOL_FOURIER 1E-8
#define KMAX  10
#define PMAX  10
#define MYPI  3.141592653589793238462643

#endif /* msm4g_constants_h */
