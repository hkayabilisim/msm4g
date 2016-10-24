/** @file msm4g_lib.h
 * @brief The declarations of the functions.
 */
#ifndef MSM4G_LIB_H
#define MSM4G_LIB_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>
#include "msm4g_types.h"

/** @brief Allocates a new dense grid.
 *
 * It create a new DenseGrid structure to hold all the values 
 * on the grid nodes. The allocated structure should be destroyed
 * using msm4g_grid_dense_destroy function.
 *
 * @param[in] h   The lattice spacing.
 * @param[in] nx  The number of nodes along x-axis.
 * @param[in] ny  The number of nodes along y-axis.
 * @param[in] nz  The number of nodes along z-axis.
 *
 * @return A new DenseGrid structure.
 */
DenseGrid *msm4g_grid_dense_new(double h, int nx, int ny, int nz);

/** @brief Reset all values to zero.
 *
 * @param[in,out] grid A grid to be reset.
 */
void msm4g_grid_dense_reset(DenseGrid *grid);

/** @brief Deallocates a given DenseGrid structure.
 *
 * Every DenseGrid structure should be deallocated by using this function to
 * properly deallocate the internal structures so that there is no memory leak.
 *
 * @param[in,out] densegrid A pointer to the DenseGrid strucuture to be deallocated.
 */
void msm4g_grid_dense_destroy(DenseGrid **densegrid);

/** @brief C1 smoothing function.
 
 \f[
 \gamma(\rho) = \left\{
 \begin{array}{ll}
 (3- \rho^2)/2, & \rho \le 1 \\
 1/\rho, & \rho \ge 1.
 \end{array}
 \right.
 \f]
 
 * @b Example:
 * 
 * To calculate \f$ \gamma(1) \f$ and \f$ \gamma'(1) \f$:
 * 
 * @code
 * double y      = msm4g_smoothing_C1(1.0, 0);
 * double yprime = msm4g_smoothing_C1(1.0, 1);
 * @endcode
 *
 * @param[in] rho        Input variable \f$ \rho \f$
 * @param[in] derivative Use the derivative of specified order.
 *
 * @return Returns \f$ \gamma(\rho) \f$.
 *
 */
double msm4g_smoothing_C1(double rho,int derivative);

/** @brief C2 smoothing function.
 
 \f[
 \gamma(\rho) = \left\{
 \begin{array}{ll}
 (15-10\rho^2 +3\rho^4)/8, & \rho \le 1 \\
 1/\rho, & \rho \ge 1.
 \end{array}
 \right.
 \f]
 *
 * @b Example:
 *
 * To calculate \f$ \gamma(1) \f$ and \f$ \gamma'(1) \f$:
 *
 * @code
 * double y      = msm4g_smoothing_C2(1.0, 0);
 * double yprime = msm4g_smoothing_C2(1.0, 1);
 * @endcode
 *
 * @param[in] rho        Input variable \f$ \rho \f$
 * @param[in] derivative Use the derivative of specified order.
 *
 * @return Returns \f$ \gamma(\rho) \f$.
 */
double msm4g_smoothing_C2(double rho,int derivative);

/** @brief C3 smoothing function.
 
 \f[
 \gamma(\rho) = \left\{
 \begin{array}{ll}
 (35- 35\rho^2 +21\rho^4 - 5\rho^6)/16, & \rho \le 1 \\
 1/\rho, & \rho \ge 1.
 \end{array}
 \right.
 \f]
 *
 * @b Example:
 *
 * To calculate \f$ \gamma(1) \f$ and \f$ \gamma'(1) \f$:
 *
 * @code
 * double y      = msm4g_smoothing_C3(1.0, 0);
 * double yprime = msm4g_smoothing_C3(1.0, 1);
 * @endcode
 *
 * @param[in] rho        Input variable \f$ \rho \f$
 * @param[in] derivative Use the derivative of specified order.
 *
 * @return Returns \f$ \gamma(\rho) \f$.
 */
double msm4g_smoothing_C3(double rho,int derivative);

/** @brief Short-range force and potential energy calculation.
 * 
 * - For each bin in the linked list
 *   - Calculate the short-range forces between the particles in the bin
 *   - Calculate the pairwise force calculatin between the neighbor bins.
 *
 * @param[in,out] binlist            The list of bins.
 * @param[in]     threshold          The range of the short-range force.
 * @param[in]     smoothing_function Pointer to the smoothing function.
 *
 * @return Total short-range potential energy in the list of bins.
 */
double msm4g_force_short(LinkedList *binlist,double threshold,msm4g_smoothing_handler smoothing_function);

/** @brief Calculates short-range forces and potential energy.
 *
 * Given a list of particles, a cut-off parameter and smoothing function
 * it calculates the short-range force interactions within the particles.
 *
 * @param[in,out] particles          The list of particles.
 * @param[in]     threshold          Cut-off parameter.
 * @param[in]     smoothing_function The smoothing function handler.
 *
 * @return Short-range potential energy in the Bin.
 */
double msm4g_force_short_withinBin(LinkedList *particles, double threshold,msm4g_smoothing_handler smoothing_function);

/** @brief Calculates short-range interactions between two set of particles.
 *
 * Given two list of particles, a cut-off parameter and smoothing function
 * it calculates the short-range interactions between two set of particles.
 *
 * @param[in,out] particlesI         A list of particles.
 * @param[in,out] particlesJ         Another list of particles.
 * @param[in]     threshold          Cut-off parameter.
 * @param[in]     smoothing_function The smoothing function handler.
 *
 * @return Short-range potential energy between the bins.
 */
double msm4g_force_short_betweenBin(LinkedList *particlesI, LinkedList *particlesJ, double threshold,msm4g_smoothing_handler smoothing_function);

/** @brief Short-range force and potential energy for a pair of particles.
 *
 * @param[in,out] particleI          The first particle.
 * @param[in,out] particleJ          The second particle.
 * @param[in]     threshold          Cut-off parameter.
 * @param[in]     smoothing_function The smoothing function handler.
 * 
 * @return Short-range potential energy.
 */
double msm4g_force_short_particlePair(Particle *particleI, Particle *particleJ, double threshold, msm4g_smoothing_handler smoothing_function);

/** @brief Sets the elements of a 3-element vector.
 *
 * The elements of a 3DVector is set to the given
 * values of x,y, and z.
 *
 * @param[in,out] d3vector A pointer to a 3-element double vector.
 * @param[in]     x        The first element.
 * @param[in]     y        The second element.
 * @param[in]     z        The third element.
 */
void msm4g_d3vector_set(D3Vector *d3vector,double x,double y,double z);

/** @brief Implemention of z = x + a*y.
 *
 * A very simple DAXPY (BLAS) implementation. The first argument z,
 * is the output. So it should be allocated before calling the function.
 *
 * @param[in,out] z The z vector in z = x + a*y.
 * @param[in]     x The x vector.
 * @param[in]     a The scalar in front of y vector.
 * @param[in]     y The y vector.
 */
void msm4g_d3vector_daxpy(D3Vector *z,D3Vector *x,double a,D3Vector *y);

/** @brief The Euclidean norm of a 3-dimensional double vector.
 *
 * This function simple calculates the Euclidean norm a vector:
 \f[
 ||x||=\sqrt{\sum_{i=1}^3 x^2_i}
 \f]
 *
 * @b Example:
 *
 * @code
 * D3Vector x;
 * double norm;
 *
 * msm4g_d3vector_set(&x,3.0,0.0,4.0);
 * norm = msm4g_d3vector_norm(&x)
 * // norm variable is 5.0 now.
 * @endcode
 *
 * @param[in] x A 3-dimensional double vector.
 * 
 * @return The norm of the vector.
 */
double msm4g_d3vector_norm(D3Vector *x);

/** @brief Calculates norm square of a vector.
 *
 * @param[in] x 3-element double vector x.
 * @return Returns \f$ ||x||^2 \f$
 *
 * @b Example:
 *
 * @code
 * D3Vector x;
 * double normsquare;
 *
 * msm4g_d3vector_set(&x,3.0,0.0,4.0);
 * normsquare = msm4g_d3vector_normsquare(&x)
 * // norm variable is 25.0 now.
 * @endcode
 */
double msm4g_d3vector_normsquare(D3Vector *x);

/** @brief Print a 3D vector.
 *
 * @param[in] x A 3-dimensional double vector.
 */
void msm4g_d3vector_print(D3Vector *x);

/** @brief Sets the elements of a 3-element vector.
 *
 * The elements of a 3D-integer vector is set to the given
 * values of x,y, and z.
 *
 * @param[in,out] i3vector A pointer to a 3-element integer vector.
 * @param[in]     x        The first element.
 * @param[in]     y        The second element.
 * @param[in]     z        The third element.
 */
void msm4g_i3vector_set(I3Vector *i3vector,int x,int y,int z);

/** @brief Copy the content of a I3Vector to another.
 *
 * @param[in,out] to    The updated vector.
 * @param[in]     from  The source vector.
 */
void msm4g_i3vector_copy(I3Vector *to, I3Vector from);

/** @brief Compare two 3D integer vector.
 *
 * @param[in] x The first vector.
 * @param[in] y The second vector.
 *
 * @return True if the vectors are same, false otherwise.
 */
Boolean msm4g_i3vector_isequal(I3Vector *x,I3Vector *y);

/** @brief Print a 3D int vector in [%d,%d,%d] format
 *
 * @param[in] x A 3-element integer array
 */
void msm4g_i3vector_print(I3Vector *x);

/** @brief Creates an empty linked list.
 *
 * This function allocates an empty linked list and returns
 * its pointer to the caller. A linked list created with this
 * function should be deallocated at the end of the program
 * by calling `msm4g_linkedlist_destroy` function.
 *
 * In order to add data to the newly creaded linked list,
 * one may use `msm4g_linkedlist_add` function.
 *
 * @return The pointer of the new linked list.
 */
LinkedList *msm4g_linkedlist_new();

/** @brief Adds a data pointer to a linked list.
 *
 * This function is used to add a data pointer to the end
 * of a linked list. Please note that, only the pointer to
 * the data is copied into the linked list. The actual data content
 * should be prepared before calling this function. And for
 * the same reason, when the linked list is deallocated, the original
 * data is not touched. Deallocation of the data should be handled
 * separately.
 *
 * @param[in,out] list The linked list to be updated.
 * @param[in]     data Pointer to the data.
 */
void msm4g_linkedlist_add(LinkedList *list,void *data);

/** @brief Get a specific element from a linked list.
 * 
 * In some situations, one may require to get a specific element
 * from the linked list by giving its exact location in the list.
 *
 * @param[in] list  The linked listed to be searched.
 * @param[in] index The zero-based index of the element.
 *
 * @return The void pointer of the requested data. 
 * If there is no data, then NULL pointer is returned.
 */
void *msm4g_linkedlist_get(LinkedList *list, int index);

/** @brief Determines the number of elements in a linked list.
 *
 * By traversing from head element to the tail, this function
 * counts the number of elements of a given linked list.
 *
 * @param[in] list The linked list whose size is to be determined.
 *
 * @return The number of elements inside the linked list.
 */
int msm4g_linkedlist_size(LinkedList *list);

/** @brief Deaallocate the linked list.
 * 
 * Since the elements of a linked list are dynamically allocated,
 * the linked list should be destroyed if it is not going to be used
 * anymore. This function deallocates all the elements (`LinkedListElement')
 * in the list.
 *
 * @warning This function does not deallocate the actual data pointed 
 * in the elements. This is the responsibility of the developer to track
 * down the location of the actual data and deallocate it at the end of the
 * program. 
 *
 * @param[in,out] list The list to be deallocated.
 */
void msm4g_linkedlist_destroy(LinkedList *list);

/** @brief Destroy the linked list including the data.
 *
 * Almost the same with the function msm4g_linkedlist_destroy
 * except that the data contained in the linked list elements 
 * are also deallocated.
 *
 * @param[in,out] list The list to be destroyed.
 */
void msm4g_linkedlist_destroyWithData(LinkedList *list);

/** @brief Allocates a brand-new empty Particle.
 *
 * This functions is the only place inside which a memory allocation
 * is carried to create a new Particle. This function holds a static Particle index
 * to keep track the newly born particles.
 *
 * @return A pointer to the newly allocated Particle structure.
 */
Particle *msm4g_particle_empty();

/** @brief Creates a new Particle with given properties.
 *
 * @param[in] mass     The mass of the Particle.
 * @param[in] location 3-element double array of location.
 * @param[in] velocity 3-element double array of velocity.
 *
 * @return The pointer of the newly allocated ::Particle.
 */
Particle *msm4g_particle_new(double mass,double *location,double *velocity);

/** @brief Creates an array of andom particles.
 *
 * This function creates an array of `Particle` structures whose properties
 * are randomly created in the range [0-1]. Therefore, the particles reside
 * in three-dimensional unit cube.
 *
 * @return The pointer to the first element of the allocated array.
 */
Particle **msm4g_particle_rand(int n);

/** @brief Load the Particle configuration from a text file.
 *
 * @param[in] filename The name of the file containing the Particle properties.
 *
 * @return The list of particles in a linked list container.
 */
LinkedList *msm4g_particle_read(const char *filename);

/** @brief Sets all of the properties of a Particle to zero value.
 *
 * Sometimes, it is desirable to use already created particles, but with
 * different properties. Then this function becomes handy, to reset its
 * properties to zero value.
 *
 * @param[in,out] Particle The Particle whose properties to be reset.
 *
 * @return The pointer of the input Particle.
 */
Particle *msm4g_particle_reset(Particle *Particle);

/** @brief Print the properties of a Particle.
 *
 * The mass, location, velocity, and other properties if it has any
 * are printed to the standard output.
 *
 * @param[in] Particle The Particle whose properties are to be printed.
 */
void msm4g_particle_print(Particle *Particle);

/** @brief Print all particles in a linked list.
 * 
 * This functions assumes that the linked list 
 * elements point to Particle data types.
 *
 * @param[in] Particlelist The list of particles to be printed.
 */
void msm4g_particle_printlist(LinkedList *Particlelist);

/** @brief Deallocates a single Particle object
 *
 * @param[in,out] Particle The Particle to be deallocated.
 */
void msm4g_particle_destroy(Particle *Particle);

/** @brief Deallocates an array of Particle objects.
 *
 * @param[in,out] Particlearray The array to be deallocates.
 * @param[in]     length    The length of the array.
 */
void msm4g_particle_destroyarray(Particle **Particlearray,int length);

/** @brief Create a unit cube for simulation 
 * 
 * Allocates a new simulation box whose location is at the origin
 *
 * @warning It is the responsibility of the caller to deallocate the
 * box by using msm4g_box_destory function.
 *
 * @return The new simulation box.
 */
SimulationBox *msm4g_box_new();

/** @brief Update the simulation box for a given list of particles.
 *
 * It creates a rectangular there dimensional box to include
 * all of the particles given in the argument. The width of the box
 * is also enlarged `margin` percent in each dimension to give some
 * slack to the particles so that they will stay in the box after time
 * integration.
 *
 * @todo If a Particle escapes from the box, the box should be updated.
 * But there is a question of how to keep the interior grid points
 * intact. As soon as the box is enlarged, new grid points should be
 * created. MSM may start from scratch whenever the box is updated,
 * but it seems quite waste of time.
 *
 * @param[in,out] box    The rectangular simulation box.
 * @param[in]     list   The list of particles.
 * @param[in]     margin The box is enlarged `margin` percent (0.10 default).
 */
void msm4g_box_update(SimulationBox *box,LinkedList *list,double margin);

/** @brief Translate the box and the particles by `delta` vector.
 *
 * @param[in,out] box    The simulation box.
 * @param[in,out] particles The list of particles.
 * @param[in]     delta  Translation amoun.
 */
void msm4g_box_translate(SimulationBox *box, LinkedList *particles,D3Vector delta);

/** @brief Translate the box and the particles so that the location of box is zero vector.
 *
 * This function is equivalant to msm4g_box_translate with delta vector is opposite to the
 * location of the box.
 *
 * @param[in,out] box    The simulation box.
 * @param[in,out] particles The list of particles.
 */
void msm4g_box_translateToOrigin(SimulationBox *box, LinkedList *particles);

/** @brief Print the details of the simulation box.
 *
 * The location, width and other features of the simulation
 * box is printed to the standard output.
 *
 * @param[in] box The box whose information is to be printed.
 */
void msm4g_box_print(SimulationBox *box);

/** @brief Deallocates the simulation box
 * 
 * The simulation box can be deallocated by using this function.
 *
 * @param[in] box The simulation box to be deallocated.
 */
void msm4g_box_destroy(SimulationBox *box);

/** @brief Create a new Bin for a given index.
 *
 * @param[in] index The new bin will have this index.
 */
Bin *msm4g_bin_new(I3Vector index);

/** @brief Generate bins for a given simulation box and Particle list.
 *
 * The particles are supposed to be contained in the box. And also the
 * the location of the simulation box is expected to be [0,0,0]. If not
 * the caller should translate the box and the particles so that the location 
 * of the box is at the origin.
 *
 * @param[in] box       The simulation box.
 * @param[in] particles The list of particles.
 * @param[in] width     The width of the bins.
 * @return The list of allocated bins.
 */
LinkedList *msm4g_bin_generate(SimulationBox *box,LinkedList *particles,double width);

/** @brief Find the neigbhors of each Bin.
 * 
 * - For each Bin in the linked list.
 *   - Check if there is a Bin in one of the 26 neighbor cubes.
 *   - If there is a bin, then add it to the neighbor list.
 *
 * @param[in,out] binlist The linked list of Bin's.
 */
void msm4g_bin_findneighbors(LinkedList *binlist);

/** @brief Returns the bin with given index vector.
 *
 * This function loops over binlist to find a bin with same
 * index information given in the argument.
 *
 * @param[in] binlist The list of bins.
 * @param[in] index   The 3-dimensional index vector.
 *
 * @return The pointer to the selected bin. Returns NULL if there is no such bin.
 */
Bin *msm4g_bin_searchByIndex(LinkedList *binlist,I3Vector index);

/** @brief Print the contents of a bin.
 * 
 * @param[in] bin The bin to be printed.
 */
void msm4g_bin_print(Bin *bin);

/** @brief Print all bins in a linked list.
 * 
 * @param[in] binlist The list of bin's.
 */
void msm4g_bin_printlist(LinkedList *binlist);

/** @brief Deallocates a bin list.
 *
 * - For each Bin in the list
 *    - Deallocate the neighbor list (without data)
 *    - Deallocate the Particle list (without data)
 * - Deallocate the bin list
 *
 * @param[in,out] binlist The linked list of Bin's.
 */
void msm4g_bin_destroy(LinkedList *binlist);
/*
 * Force calculation
 */
void msm4g_acceleration(double *a,double *r,int n,int d,double *m,double G);

/*
 *Energy calculations
 */
void msm4g_energy(double *pot,double *kin,double *tot,double *r,double *v,int n,int d,double *m,double G);

/*
 * Time integration schemes
 */
void msm4g_forwardeuler(double *r,double *v,double *a,double dt,int n,int d,double *m,double G);

/** @brief Leap-Frog time-integration scheme.
 *
 * @param[in,out] r  The locations of the particles.
 * @param[in,out] v  The velocities of the particles.
 * @param[in,out] a  The acceleration of the particles.
 * @param[in,out] a1 The acceleration of the particles on half-step time.
 * @param[in]     dt The timestep.
 * @param[in]     n  The number of particles.
 * @param[in]     d  The spatial dimension.
 * @param[in]     m  The mass of the particles.
 * @param[in]     G  the gravitational constant.
 */
void msm4g_leapfrog(double *r,double *v,double *a,double *a1,double dt,int n,int d,double *m,double G);
/*
 * IO-related
 */
void msm4g_print_energy(double pot0,double kin0,double tot0,double pot,double kin,double tot);

/*
 * Memory related
 */
void msm4g_zeros(double *x,int n);

/** @brief n-degree generalized Cantor pairing function.
 *
 * This method is the direct implementation of 
 * the n-degree generalized Cantor pairing function
 * defined in Theorem 2.1 in Lisi, M. (2007). Some remarks on 
 * the Cantor pairing function. Le Matematiche, 62(1), 55–65.
 *
 * \f[
  <x_1,\ldots,x_2> = \sum_{h=1}^n\left\{ \frac{1}{h!}\prod_{j=0}^{h-1} \left[ \left( \sum_{i=1}^hx_i\right)+j\right] \right\}
 \f]
 *
 * @param[in] x An integer vector in N^n.
 * @param[in] n The dimension.
 *
 * @return The output of n-degree generalized Cantor pairing function.
 */
int msm4g_math_cantor(int *x,int n);

/** @brief Factorial function.
 *
 * @param[in] n A natural number.
 * @return factorial of n.
 */
int msm4g_math_factorial(int n);

#endif /* MSM4G_LIB_H */

