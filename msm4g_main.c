/** @file msm4g_main.c
 * @brief The main entry point of the MSM4G application.
 *
 * @mainpage
 *
 * MSM4G program is aimed to use the Multi-Level Summation Method (MSM) for (4)
 * Gravitational (G) N-Body problems @anchor logonote *.
 *
 * @section datastructures Data Structures
 *
 * - Particle is used to represent a celestial object.
 * - LinkedList is a generic two-way linked listed implementation.
 * - Bin is a container to hold the particles in geometric hashing.
 * - SimulationBox is a rectangular box enclosing all of the particles.
 *
 * @section shortrange Short-Range Calculations
 * 
 * The simulation domain is divided into cubic boxes called Bin whose edges are
 * set to the cut-off parameter. Depending on the particle distribution and the 
 * arrangement of the simulation box, there might be considerable amount of Bin's
 * containing no particles. In order to save some space, the Bin's are stored as
 * a linked list. Only the Bin's containing at least one particle are stored in the
 * linked list.
 * 
 * The neighbors of a specific Bin are also stored in a linked list. Again, a Bin
 * is added as a neighbor only if it contains at least a particle.
 *
 * In order to exploit the Newton's third law during the short-range calculations,
 * each Bin is enumerated by a unique integer corresponding to its index in the
 * three dimensional simulation box. An exotic choise is Cantor's generalized
 * pairing function (See. msm4g_math_cantor()). Whenever a new Bin is allocated, its 
 * cantor index is calculated and stored in the Bin structure. Since the Bin structure
 * already has the index information, there is no need to invert the cantor index.
 *
 * @section tests Test Functions
 * 
 * A user who wants to use MSM4G in his/her application or
 * tries to understand the details of the package may prefer
 * to look at the test functions in msm4g_tests.c as they are
 * self-contained and test various features individually.
 *
 * @ref logonote "*" The artwork at top of the page is the picture of the Andromeda Galaxy
 * which is taken using amateur refractor telescope. Courtesy of Kanwar Singh. 
 * The link to the higher resolutions is: https://commons.wikimedia.org/wiki/File:ANDROMEDA_GALAXY.jpg
 */

#include "msm4g_lib.h"
#include "msm4g_tests.h"


int main()
{
    msm4g_unit_test_all();
    
    return 0;
}
