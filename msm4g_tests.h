/** @file msm4g_tests.h
 * @brief The declarations of the test functions.
 *
 * The test functions are self-contained so they don't require
 * inputs. They only return a boolean indicating the fate of the
 * test. The test function msm4g_unit_test_all() runs all of the tests.
 */
#ifndef MSM4G_TESTS_H
#define MSM4G_TESTS_H

#include "msm4g_lib.h"

/** @brief Runs all tests.
 *
 * All test functions are sequentiall run and their results
 * are printed to standard output. If one of the tests is
 * failed, then the corresponding test information is printed
 * to the standard error.
 */
void msm4g_unit_test_all();

/** @brief Checks the size() functionality of linked lists.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_1();

/** @brief Traversing the linked list from tail to head.
 *
 * Using the `prev` pointer element in the `LinkedListElement`
 * one can travel from the tail element to the head.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_2();

/** @brief Traversing the list from head to tail.
 *
 *  @return true if the test is succesfull, false otherwise.
 */

Boolean msm4g_unit_test_3();
/** @brief Reading particle configurations from a text file.
 *
 * The masses, locations and velocities of three particles are read from
 * eight.ini file and compared to the true values. If there is a slight
 * variation, the test returns with fail.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_4();

/** @brief Creating a box surrounding the particles.
 *
 * Given a collection of particles, this test determines a rectangular
 * simulation box surrounding all of them and compares the result
 * with the expected box. In this example, 10 random particles are created
 * in unit 3D cube. Then the location of one of them is set to [-1,-2,-3].
 * Another one is set to [1,2,3]. Given this configuration the location,
 * and the width of the box shoould be [-1,-2,-3] and [2,4,6] respectively.
 * But then a 0.5 margin is also used to enlarge the box. Then the new width
 * should be [3,6,9]. The location should also be shifted so that the there is
 * a margin on both minus and plus axis. So the new location is [-1.5,-3,-4.5].
 * At the end of the test, the calculated location and width is compared with the
 * expecteds.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_5();

/** @brief Generating bins for the particles in eight.ini file.
 *
 * The flow of this test is as follows:
 * 
 * - Loading the particles from eight.ini file
 * - Creating a simulation box around the particles with 0.5 extra margin
 * - Generating the bins
 *
 * @verbatim
 30  +----------+----------+----------+     Bin Width: 10
 |   |          |  4       |          |
 |   |          |   5      |          |     Particle Location
 |   |          |    6     |          |     -------- --------
 |   |          |          |          |     1        [15,15, 5]
 20  +----------+----------+----------+     2        [15,12, 5]
 |   |          |          |          |     3        [29, 2,25]
 |   |          |     1    |          |     4        [12,26, 5]
 |   |          |          |          |     5        [13,25, 5]
 |   |          |     2    |          |     6        [14,24, 5]
 10  +----------+----------+----------+
 |   |          |          |          |     Bins
 |   |          |          |          |     ----
 |   |          |          |          |     i j k   Particles Neighbor Bins
 |   |          |          |         3|     - - -   --------- -------------
 0   +----------+----------+----------+     1 1 0   1,2       [1,2,0]
                                            2 0 2   3         none
     0----------10---------20---------30    1 2 0   4,5,6     [1,1,0]
        i=0         i=1         i=2

  @endverbatim
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_6();

/** @brief Testing Cantor pairing function.
 *
 * The test function is used to check the output of 
 * Cantor pairing function msm4g_math_cantor() applied to a 4x4 matrix.
 * If the enumeration is not equal to the one below, 
 * the test returns false.
 * The test also checks the output of msm4g_math_factorial() function 
 * for the integers, 0, 1, 2, 3,  4, and 5.
 *
 * i/j | j=0 | j=1 | j=2 | j=4 |
 * --- | --- | --- | --- | --- |
 * i=0 |  0  |  1  |  3  |  6  |
 * i=1 |  2  |  4  |  7  |  11 |
 * i=2 |  5  |  8  |  12 |  17 |
 * i=3 |  9  |  13 |  18 |  24 |
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_7();

/** @brief Testing smoothing functions.
 *
 * Checks the implementation of the smoothing functions and their
 * derivatives on p=0.0 and p=1.0;
 * The \f$ \gamma \f$ functions tested are msm4g_smoothing_C1(),
 * msm4g_smoothing_C2(), and msm4g_smoothing_C3().
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_8();

/** @brief Short-range energy calculation test.
 *
 * The setup of the particles and the simulation box is the same
 * with msm4g_unit_test_6(). This time short-range force and energy calculation
 * is carried out.
 * 
 * The manually calculated short-range potential regarding to the setup
 * is 1.51303181237454868855 which is compared to the calculated potential energy.
 * If the relative error in the potential energy is greater than 1E-15, the test 
 * is failed.
 *
 * @note In this test, C1 gamma smoothing function msm4g_smoothing_C1() is used with
 * cut-off parameter a=10.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_9();

/** @brief Testing DenseGrid implementation.
 * 
 * This function tests the basic functionality of DenseGrid related 
 * functions.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_10();

/** @brief Testing anterpolation.
 *
 * In this test, 8 unit mass particles are created on the following locations.
 * The degree of base polynomials is set to p=3. Then a simulation box
 * is created with margin=0. Since p=3, an extra padding is added around the
 * boundary as shown in the following figure. Then the anterpolation is run.
 * Since Hermite nodal basis functions are used in the anterpolation, the grid
 * masses on the 8 interior nodes will be 1, and the rest should be zero.
 * 
 @verbatim
 30  +----------+----------+----------+     Lattice spacing (h): 10
 |   |          |          |          |
 |   |          |          |          |     Particle Location
 |   |          |          |          |     -------- --------
 |   |          |          |          |     1        [10,10,10]
 20  +----------O----------O----------+     2        [10,10,20]
 |   |          |          |          |     3        [10,20,10]
 |   |          |          |          |     4        [10,20,20]
 |   |          |          |          |     5        [20,10,10]
 |   |          |          |          |     6        [20,10,20]
 10  +----------O----------O----------+     7        [20,20,10]
 |   |          |          |          |     8        [20,20,20]
 |   |          |          |          |
 |   |          |          |          |
 |   |          |          |          |
 0   +----------+----------+----------+
 
     0----------10---------20---------30
 @endverbatim
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_11();

/** @brief To be or not to be!
 * 

 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_12();

#endif /* MSM4G_TESTS_H */
