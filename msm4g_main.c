/** @file msm4g_main.c
 * @brief The main entry point of the MSM4G application.
 *
 * A user who wants to use MSM4G in his/her application or 
 * tries to understand the details of the package may prefer
 * to look at the test functions in msm4g_tests.c as they are
 * self-contained and test various features individually.
 *
 * @mainpage
 *
 * MSM4G program is aimed to apply the Multi-Level Summation Method (MSM) for (4)
 * Gravitational (G) N-Body problems. 
 *
 * @section datastructures Data Structures
 *
 * - Particle is used to represent a celestial object.
 * - LinkedList is a generic two-way linked listed implementation.
 * - Bin is a container to hold the particles in geometric hashing.
 * - SimulationBox is a rectangular box enclosing all of the particles.
 *
 */

#include "msm4g_lib.h"
#include "msm4g_tests.h"


int main()
{
    msm4g_unit_test_all();
    
    return 0;
}
