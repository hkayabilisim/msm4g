/** @mainpage
 *
 * @section Introduction
 *
 * MSM4G program is aimed to apply the Multi-Level Summation Method (MSM) for (4)
 * Gravitational (G) N-Body problems. 
 *
 * @section datastructures Data Structures
 *
 * - Body is used to represent a celestial object.
 * - LinkedList is a generic two-way linked listed implementation.
 * - Bin is a container to hold the bodies in geometric hashing.
 * - SimulationBox is a rectangular box enclosing all of the bodies.
 *
 */

#include "msm4g_lib.h"
#include "msm4g_tests.h"

int main()
{
    msm4g_unit_test_all();
    
    return 0;
}
