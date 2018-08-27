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
#include <assert.h>
#include "msm4g_types.h"
#include "msm4g_bases.h"
#include "msm4g_constants.h"

#define MSM4G_MAX3(A,B,C) (((A)>(B)) ? (((A)>(C)) ? (A) : (C)) : (((B)>(C)) ? (B) : (C)) )

/** @brief Calculates the stencils
 *
 * This function calculates the stencil values
 * corresponding to the specified level. It modifies
 * the stencil located in the simulation object.
 *
 * @param[in] simulation Simulation object
 * @param[in] level      Level in the range [1,L+1]. l = 1 corresponds
 * to the lowest level stencil, whereas l = L + 1 is for Fourier
 * stencil.
 */
void msm4g_stencil(Simulation *simulation, int level);

/** @brief Creates a new simulation object
 * 
 * It reads particles from a file and builds a
 * simulation object.
 *
 * @param[in] datafile Filename containing the particles
 * @param[in] box      Simulation box object
 * @param[in] periodic If true, periodics boundary condition applies
 * @param[in] order    Order of accuracy. It can only be 4 or 6
 * @param[in] abar     Relative-cutoff
 * @param[in] mu       Quasi-interpolation parameter
 * @todo Change "4 or 6" statement when it is needed.
 * @return Simulation object
 */
Simulation *msm4g_simulation_new(char *datafile,SimulationBox *box,Boolean periodic,int order,double abar,int mu);

/** @brief Executes a simulation
 *
 * @param[in] simulation Pointer to the simulation object
 */
void msm4g_simulation_run(Simulation *simulation);

/** @brief Deallocates a given simulation object
 *
 * @param[in] simulation Pointer to the simulation object
 */
void msm4g_simulation_delete(Simulation *simulation);

/** @brief Anterpolation
 * @image html box.png
 */
void msm4g_anterpolation(Simulation *simulation);

/** @brief Restriction
 * It performs restriction step from level l-1
 * to l  where l is in the range [2, L].
 * @param[in] simulation Simulation object
 * @param[in] l          Level number in [2,L]
 */
void msm4g_restriction(Simulation *simulation, int l);

/** @brief Prolongation
 * Execute all prolongation steps.
 * @param[in] simulation Simulation object
 */
void msm4g_prolongation(Simulation *simulation);

/** @brief Interpolation 
 * This function calculates the long-range acceralation of the
 * particles by interpolating the grid potentials on the finest level and
 * the B-Splines.
 *
 * @param[in] simulation the object containing all information about the simulation
 */
void msm4g_interpolation(Simulation *simulation);

/** @brief Print the content of grid to the stdout.
 *
 * @param[in] grid A grid to be printed.
 */
void msm4g_grid_print(AbstractGrid *grid);

/** @brief Deallocates a grid.
 *
 * This function first calls the destructor of the 
 * grid, and then deallocates itself.
 *
 * @param[in,out] grid The grid to be allocated.
 */
void msm4g_grid_destroy(AbstractGrid **grid);

/** @brief A DenseGrid constructor.
 *
 * It create a new DenseGrid structure to hold all the values 
 * on the grid nodes. The internal data structure is later dealocated
 * in the msm4g_grid_dense_destroy() function.
 *
 * @param[in] nx  The number of nodes along x-axis.
 * @param[in] ny  The number of nodes along y-axis.
 * @param[in] nz  The number of nodes along z-axis.
 * @param[in] hx  The lattice spacing along x-axis.
 * @param[in] hy  The lattice spacing along y-axis.
 * @param[in] hz  The lattice spacing along z-axis.
 *
 * @return A new DenseGrid structure.
 */
AbstractGrid *msm4g_grid_dense_new(int nx, int ny, int nz,double hx, double hy, double hz);

/** @brief Sets a value to a specific coordinate of the DenseGrid.
 *
 * With this function, one can change a single element of DenseGrid.
 *
 * @param[in,out] grid  A DenseGrid whose element is changed.
 * @param[in]     i     The coordinate in the x-axis.
 * @param[in]     j     The coordinate in the y-axis.
 * @param[in]     k     The coordinate in the z-axis.
 * @param[in]     value The new value of the element.
 */
void msm4g_grid_dense_setElement(AbstractGrid *grid,int i,int j,int k,double value);


/** @brief Gets a single element from a DenseGrid object.
 *
 * Returns the (i,j,k) element of the dense matrix.
 *
 * @param[in] grid The AbstractGrid object.
 * @param[in] i    The coordinate in the x-axis.
 * @param[in] j    The coordinate in the y-axis.
 * @param[in] k    The coordinate in the z-axis.
 *
 * @return The value of the element corresponding to (i,j,k) coordinate.
 */
double msm4g_grid_dense_getElement(AbstractGrid *grid,int i,int j,int k);

/** @brief Calculates inner product of two dense grids.
 *
 * Inner product of two dense grid is the sum of pairwise
 * multiplications of the dense grids.
 *
 * @param[in] grid1 First dense grid
 * @param[in] grid2 Second dense grid
 * @return the inner product of the grids.
 */
double msm4g_grid_dense_innerProduct(AbstractGrid *grid1,AbstractGrid *grid2);

/** @brief Sum of grid elements
 * @param[in] grid A dense grid
 * @return sum of grid elements
 */
double msm4g_grid_dense_sum(AbstractGrid *grid);

/** @brief Add two dense grids
 * It implements the following sum:
 * grid1 = grid1 + grid2
 * @param[in,out] grid1 First dense grid
 * @param[in] grid2 Second dense grid
 */
void msm4g_grid_dense_add(AbstractGrid *grid1,AbstractGrid *grid2);

/** @brief Reset all values in the grid to the given value.
 *
 * This function resets the elements of the internal
 * dense matrix to the given value.
 *
 * @param[in,out] grid   A grid to be reset.
 * @param[in]     value  The grid is to be reset to value.
 */
void msm4g_grid_dense_reset(AbstractGrid *grid,double value);

/** @brief Deconstructor for the DenseGrid object.
 *
 * Internal data structure used in a DenseGrid object is deallocated by using 
 * this function. The deallocation of the grid object itself is handled in 
 * msm4g_grid_destroy() function.
 *
 * @param[in,out] densegrid A pointer to the DenseGrid strucuture to be deallocated.
 */
void msm4g_grid_dense_destroy(AbstractGrid **densegrid);

/** @brief Calculate grid-to-grid interaction to get grid potentials
 *
 * @todo Write down the LateX formula
 *
 * @param[in] stencil
 * @param[in] gridmass
 * @param[in,out] gridpotential
 */
void msm4g_grid_potential(AbstractGrid *stencil, AbstractGrid *gridmass, AbstractGrid *gridpotential);

/** @brief Even-powered Softener
 *
 \f[
 \gamma(\rho) = \left\{
 \begin{array}{ll}
 \tau_\nu(\rho^2), & 0 \le \rho \le 1, \\
 1/\rho, & \rho \ge 1.
 \end{array}
 \right.
 \f]
 where
 \f[
 \tau_\nu(s) = \sum_{k=0}^{\nu-1}{-1/2\choose k}(s-1)^k
 \f]
 * @param[in] rho Input value
 * @param[in] nu  order of accuracy
 *
 * @return Returns \f$ \gamma(\rho) \f$
 */
double msm4g_smoothing_gama(double rho,int nu);

/** @brief Derivative of even-powered softener
 *
 * @param[in] rho Input value
 * @param[in] nu  order of accuracy
 *
 * @return Returns \f$ \gamma'(\rho) \f$
 */
double msm4g_smoothing_gamaprime(double rho,int nu);

/** @brief Compute splitting kernels
 *
 * @todo Include the equations in LateX format
 *
 * @param[in] l Level in the range of [0,L+1]
 * @param[in] L Number of levels.
 * @param[in] r Distance value
 * @param[in] a Absolute cutoff
 * @param[in] beta Splitting parameter
 * @param[in] nu Order
 * @return the value of the kernel function
 */
double msm4g_kernel(int l,int L, double r,double a,double beta,int nu);

/** @brief Short-range force and potential energy calculation.
 * 
 * - For each bin in the linked list
 *   - Calculate the short-range forces between the particles in the bin
 *   - Calculate the pairwise force between the neighbor bins.
 *
 * @param[in,out] binlist    The list of bins.
 * @param[in]     threshold  The range of the short-range force.
 * @param[in]     simulation Simulation object.
 *
 * @return Total short-range potential energy in the list of bins.
 */
double msm4g_force_short(LinkedList *binlist,double threshold, Simulation *simulation);

/** @brief Calculates short-range forces and potential energy.
 *
 * Given a list of particles,
 * it calculates the short-range force interactions within the particles.
 *
 * @param[in,out] particles  The list of particles.
 * @param[in]     threshold  Cut-off parameter.
 * @param[in]     simulation Simulation object.
 *
 * @return Short-range potential energy in the Bin.
 */
double msm4g_force_short_withinBin(LinkedList *particles, double threshold,Simulation *simulation);

/** @brief Calculates short-range interactions between two set of particles.
 *
 * Given two list of particles,
 * it calculates the short-range interactions between two set of particles.
 *
 * @param[in,out] particlesI  A list of particles.
 * @param[in,out] particlesJ  Another list of particles.
 * @param[in]     threshold   Cut-off parameter.
 * @param[in]     simulation Simulation object.
 *
 * @return Short-range potential energy between the bins.
 */
double msm4g_force_short_betweenBin(LinkedList *particlesI, LinkedList *particlesJ, double threshold, Simulation *simulation);

/** @brief Short-range force and potential energy for a pair of particles.
 *
 * @param[in,out] particleI   The first particle.
 * @param[in,out] particleJ   The second particle.
 * @param[in]     threshold   Cut-off parameter.
 * @param[in]     simulation  Simulation object.
 * 
 * @return Short-range potential energy.
 */
double msm4g_force_short_particlePair(Particle *particleI, Particle *particleJ, double threshold, Simulation *simulation);

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
 * @param[in] x        x-axis of the location.
 * @param[in] y        y-axis of the location.
 * @param[in] z        z-axis of the location.
 *
 * @return The pointer of the newly allocated ::Particle.
 */
Particle *msm4g_particle_new(double mass,double x, double y,double z);

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
 * @param[in,out] N Number of particles
 * @return Pointer to a particle array
 */
Particle *msm4g_particle_read(const char *filename,int *N);

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
 * The mass, location, and other properties if it has any
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
 * @todo Make sure that the width is of the form 2^k * h for some k.
 *
 * @return The new simulation box.
 */
SimulationBox *msm4g_box_new();

/** @brief Creates a cube
 *
 * @param[in] location Location of the cube
 * @param[in] width    Width of the cube
 * @return Reference to the simulation box object
 */
SimulationBox *msm4g_box_newCube(double location, double width);

/** @brief Update the simulation box for a given list of particles.
 *
 * It creates a rectangular there dimensional box to include
 * all of the particles given in the argument. The width of the box
 * is also enlarged `margin` percent in each dimension to give some
 * slack to the particles so that they will stay in the box after time
 * integration.
 *
 * The width of the box is arranged to be a multiple of lattice
 * spacing h. 
 *
 * The degree of the base polynomials is also needed to make sure that
 * there is enough extra padding around the boundary.
 *
 * @todo If a Particle escapes from the box, the box should be updated.
 * But there is a question of how to keep the interior grid points
 * intact. As soon as the box is enlarged, new grid points should be
 * created. MSM may start from scratch whenever the box is updated,
 * but it seems quite waste of time.
 *
 * @param[in,out] box       The rectangular simulation box.
 * @param[in]     particles The list of particles.
 * @param[in]     margin    The box is enlarged `margin` percent (0.10 default).
 * @param[in]     h         Lattice spacing in the finest level grid.
 * @param[in]     p         The degree of the base polynomials. p should be odd and larger than 0. 
 */
void msm4g_box_update(SimulationBox *box,LinkedList *particles,double margin,double h,double p);

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
 * @param[in] n         Number of particles
 * @param[in] width     The width of the bins.
 * @return The list of allocated bins.
 */
LinkedList *msm4g_bin_generate(SimulationBox *box,Particle *particles,int n,double width);

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

/** @brief Calculates Quasi-interpolation coefficients 
 * 
 * Calculates the coefficients (referred as omegaprime in doi:10.1063/1.4943868) 
 * required in Quasi-interpolation. They depend on the order of B-Spline
 * and some adjustable parameter mu.
 * 
 * The output of this function is used to build the stencils.
 * Large mu means more accurate interpolation but at the same it increases the cost.
 *
 * The callee is responsible to deallocate the array when it is no longer needed.
 *
 * @todo Quasi-interpolation coefficients can only be determined for p=4 (Cubic) and
 * p=6 (Quintic) B-Splines. It should be written in more generic way to cover higher
 * order B-Splines.
 *
 * @param[in] mu Quasi-interpolation parameter
 * @param[in] p Order of B-Spline. It can be 4 (Cubic) and 6 (Quintic) only.
 * @return Quasi-interpolation coefficient array of size mu + p/2 + 1.
 */
double *msm4g_util_omegaprime(int mu,int p) ;

/** @brief L2-norm of a vector
 * 
 * @param[in] x input vector
 * @param[in] n size of the input vector
 * @return L2-norm of the vector x
 */
double msm4g_util_norm(double x[], int n);

/** @brief L2-norm of x-y
 *
 * @param[in] x first vector
 * @param[in] y second vector
 * @param[in] n size of the vectors
 * @return L2-norm of x-y
 */
double msm4g_util_diffnorm(double x[], double y[], int n);

/** @brief Solve Ax=y by Gaussian elimination
 *
 * This solver is very simple featuring:
 *   - Does not use any pivoting.
 *   - Overwrites the matrix A and vector y.
 *   - At the end of the process, the right-hand side vector 
 *     y becomes the solution.
 *
 * @param[in] n size of the system
 * @param[in] A matrix of size n by n
 * @param[in,out] y right hand-size vector of size n. It is also the solution.
 */
void msm4g_util_gausssolver(int n, double *A, double *y);

/** @brief Enumerate grid points on the surface of a cubic grid
 *
 * Consider a cube centered at the origin and with size of 2n x 2n x 2n.
 * The task is to enumerate the grid points on the surface of the cube.
 *
 * @param[in] n  Size of the cube is 2n x 2n x 2n
 * @param[in] sp SimulationParameters object. 
 * @return the number of grid points on the cube surface
 */
int msm4g_util_face_enumerate(int n,SimulationParameters *sp);

/** @brief Determines Ewald's splitting parameter
 *
 * Determines beta parameter such that
 * erfc(beta * aL)/aL is very small
 *
 * @todo Elaborate the scheme used in the function.
 *
 * @param[in] aL Cutoff at the coarsest grid
 * @return Ewald's splitting parameter beta
 */
double msm4g_util_choose_beta(double aL);

/** @brief Calculates the coefficients needed for interpolation at (L+1)th level
 * Refer eq:calpha in the manuscript
 * @param[in] k wavenumber in Fourier transform
 * @param[in] M The number of grid points
 * @param[in] nu Order of interpolation
 * @return value of the coefficient
 */
double msm4g_util_calculate_c(int k, double M, int nu);

/** @brief Implements the binomial function
 * @param[in] n an integer
 * @param[in] k an integer
 * @return nchoosek(n,k)
 */
double msm4g_util_nchoosek(int n,int k);

/** @brief Implements the restriction operator
 * 
 * @todo Put the formula
 *
 * @param[in] nu order of B-Splines
 * @param[in] k  an integer
 */
double msm4g_util_jn(int nu,int k);

#endif /* MSM4G_LIB_H */

