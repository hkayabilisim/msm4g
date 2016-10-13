#ifndef MSM4G_TESTS_H
#define MSM4G_TESTS_H

#include "msm4g_types.h"

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
/** @brief Reading body configurations from a text file.
 *
 * The masses, locations and velocities of three bodies are read from
 * eight.ini file and compared to the true values. If there is a slight
 * variation, the test returns with fail.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_4();

/** @brief Creating a box surrounding the bodies.
 *
 * Given a collection of bodies, this test determines a rectangular
 * simulation box surrounding all of them.
 *
 * @return true if the test is succesfull, false otherwise.
 */
Boolean msm4g_unit_test_5();

#endif /* MSM4G_TESTS_H */
