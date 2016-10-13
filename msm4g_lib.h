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


/* Body related */
Body *msm4g_newbody(double mass,double *location,double *velocity);
Body *msm4g_newzerobody();

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
void msm4g_print_body(double *r,double *v,double *a,int i);
void msm4g_print_energy(double pot0,double kin0,double tot0,double pot,double kin,double tot);
void msm4g_printbody(Body *body);

/*
 * Memory related
 */
void msm4g_zeros(double *x,int n);

#endif /* MSM4G_LIB_H */

