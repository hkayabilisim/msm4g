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

/** @brief Implemention of z = a*x + y.
 *
 * A very simple DAXPY (BLAS) implementation. The first argument z,
 * is the output. So it should be allocated before calling the function.
 *
 * @param[in,out] z The z vector in z=a*x+y.
 * @param[in]     a The scalar in front of x vector.
 * @param[in]     x The x vector.
 * @param[in]     y The z vector.
 */
void msm4g_d3vector_daxpy(D3Vector *z,double a,D3Vector *x,D3Vector *y);

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

/** @brief Creates a new body with given properties.
 *
 * @param[in] mass     The mass of the body.
 * @param[in] location 3-element double array of location.
 * @param[in] velocity 3-element double array of velocity.
 *
 * @return The pointer of the newly allocated ::Body.
 */
Body *msm4g_body_new(double mass,double *location,double *velocity);

/** @brief Creates `n` number of random bodies.
 *
 * This function creates an array of `Body` structures whose properties
 * are randomly created in the range [0-1]. Therefore, the bodies reside
 * in three-dimensional unit cube.
 *
 * @return The pointer to the first element of the allocated array.
 */
Body *msm4g_body_rand(int n);

/** @brief Sets all of the properties of a body to zero value.
 *
 * Sometimes, it is desirable to use already created bodies, but with
 * different properties. Then this function becomes handy, to reset its
 * properties to zero value.
 *
 * @param[in,out] body The body whose properties to be reset.
 *
 * @return The pointer of the input body.
 */
Body *msm4g_body_reset(Body *body);

/** @brief Print the properties of a body.
 *
 * The mass, location, velocity, and other properties if it has any
 * are printed to the standard output.
 *
 * @param[in] body The body whose properties are to be printed.
 */
void msm4g_body_print(Body *body);

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

/** @brief Update the simulation box for a given list of bodies.
 *
 * It creates a rectangular there dimensional box to include
 * all of the bodies given in the argument. The width of the box
 * is also enlarged `margin` percent in each dimension to give some
 * slack to the bodies so that they will stay in the box after time
 * integration.
 *
 * @todo If a body escapes from the box, the box should be updated.
 * But there is a question of how to keep the interior grid points
 * intact. As soon as the box is enlarged, new grid points should be
 * created. MSM may start from scratch whenever the box is updated,
 * but it seems quite waste of time.
 *
 * @param[in,out] box    The rectangular simulation box.
 * @param[in]     list   The list of bodies.
 * @param[in]     margin The box is enlarged `margin` percent (0.10 default).
 */
void msm4g_box_update(SimulationBox *box,LinkedList *list,double margin);

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
 * @param[in,out] r  The locations of the bodies.
 * @param[in,out] v  The velocities of the bodies.
 * @param[in,out] a  The acceleration of the bodies.
 * @param[in,out] a1 The acceleration of the bodies on half-step time.
 * @param[in]     dt The timestep.
 * @param[in]     n  The number of bodies.
 * @param[in]     d  The spatial dimension.
 * @param[in]     m  The mass of the bodies.
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

#endif /* MSM4G_LIB_H */

