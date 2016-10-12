/**
 *         Multi-Level Summation Method (MSM)
 *                   for
 *         Gravitational N-body Problems
 *
 *
 *             Function Declerations
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>

#include "nbody_types.h"

#define NBODY_ERR_INPUT_ISNOT_3D -1

/** @brief Running all unit tests */
void nbody_unit_test_all();
/** @brief checking the size() functionality of linked lists */
bool nbody_unit_test_1();
/** @brief checking the traversal from tail to head */
bool nbody_unit_test_2();
/** @brief checking the traversal from head to tail */
bool nbody_unit_test_3();
/** Reading body masses,locations and velocities from a text file. */
bool nbody_unit_test_4();
bool nbody_unit_test_5();


int nbody_test_eight();

void nbody_linkedlist_destroy(LinkedList *list);
void nbody_linkedlist_add(LinkedList *list,void *data);
LinkedList *nbody_linkedlist_new();
int nbody_linkedlist_size(LinkedList *list);


/* Body related */
Body *nbody_newbody(double mass,double *location,double *velocity);
Body *nbody_newzerobody();
Body *nbody_newrandbody();
Body *nbody_resetbody(Body *body);

/* Geometry */
void nbody_createBox(SimulationBox *box,LinkedList *bodies);
/*
 * Force calculation
 */
void nbody_acceleration(double *a,double *r,int n,int d,double *m,double G);

/*
 *Energy calculations
 */
void nbody_energy(double *pot,double *kin,double *tot,double *r,double *v,int n,int d,double *m,double G);
/*
 * Time integration schemes
 */
void nbody_forwardeuler(double *r,double *v,double *a,double dt,int n,int d,double *m,double G);
void nbody_leapfrog(double *r,double *v,double *a,double *a1,double dt,int n,int d,double *m,double G);
/*
 * IO-related
 */
void nbody_print_body(double *r,double *v,double *a,int i);
void nbody_print_energy(double pot0,double kin0,double tot0,double pot,double kin,double tot);
void nbody_printbody(Body *body);

/*
 * Memory related
 */
void nbody_zeros(double *x,int n);

