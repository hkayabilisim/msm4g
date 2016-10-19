/** @file msm4g_types.h
 * @brief The data type definitions required in the package.
 *
 * It is already included in @link msm4g_lib.h @endlink so that it
 * doesn't have to be explicity included by the user.
 */
#ifndef MSM4G_TYPES_H
#define MSM4G_TYPES_H

/** @brief A custom definition of boolean type.
 *
 * Although, C99 has boolean type, I prefer to stay on ANSI C.
 */
typedef enum
{
    false ,  /*!< logical false */
    true     /*!< logical true  */
} Boolean;

/** @brief 3-dimensional double vector.
 *
 * Since, the domain is 3-dimensional, it can be handy to define
 * a structure to hold three values in one entity
 */
typedef struct D3Vector
{
    double value[3];
} D3Vector;

/** @brief 3-dimensional integer vector.
 */
typedef struct I3Vector
{
    int value[3];
} I3Vector;


/** @brief An element of a LinkedList.
 *
 * A LinkedList is composed of a collection of LinkedListElement structure
 * which contains a void pointer to the data so that one can store different
 * data types in the list. In this way, the same linked list implementation and
 * its corresponding functions can be reused for all data types.
 *
 * The structure also contans `next`, and `prev` pointers to traverse in the list
 * in both directions.
 */
typedef struct LinkedListElement
{
    void *data;                     /**< Pointer to the data */
    struct LinkedListElement *next; /**< Pointer to the next element */
    struct LinkedListElement *prev; /**< Pointer to the previous element */
} LinkedListElement;

/** @brief A structure for dual linked list implementation. 
 * 
 * Our linked listed implementation contains two pointers of 
 * LinkedListElement: head and tail. These pointers only change
 * if a new element is added or removed from the list. In other cases,
 * they always point to the head and tail elements in the list. If one needs
 * to traverse the list in either direction, the developer needs to use
 * a temporary pointer so that the head and tail pointers are left intact.
 */
typedef struct LinkedList
{
    struct LinkedListElement *head;
    struct LinkedListElement *tail;
} LinkedList;

/** @brief A basic representation of a celestial body.
 *
 * Every celestial body in the gravitational n-body problem is 
 * represented with this simple Particle structure containing basic features
 * such as mass, location, velocity, and force.
 *
 * @warning The index is supposed to be read-only in the user-end.
 * This variable is created when the particle is allocated for the very first 
 * time inside the msm4g_particle_empty() function. This means that user should
 * not allocate a particle manually. He/She should always use msm4g_particle functions
 * to handle the particles.
 */
typedef struct Particle
{
    int    index;      /**< Index of the particle (see warning in the description) */
    double m;          /**< mass      */
    D3Vector r;        /**< location  */
    D3Vector v;        /**< velocity  */
    D3Vector f;        /**< force produced by short-range interactions. */
    D3Vector flong;    /**< force produced by long-range interactions.  */
} Particle;

/** @brief The cubical domains obtained by dividing the simulation box.
 * 
 * Bin is a rectangular box created by dividing the simulation box.
 * The member cantorindex is used to assign a unique number to each
 * bin. This number is calculated only during the generation of a new
 * Bin.
 */
typedef struct Bin
{
    I3Vector index;         /**< The index of the bin balong each axis */
    int      cantorindex;   /**< Cantor index corresponding to index variable */
    LinkedList *neighbors;  /**< The list of neightbor bins containing at leeast one particle */
    LinkedList *particles;  /**< The list of particles in the bin */
} Bin;

/** @brief Smoothing function handler
 *
 * A function pointer to represent smoothing functions
 * which takes a double argument and the order of
 * the derivative and returns the output
 * as a double variable.
 *
 * If the derivative is set to zero, then no derivative is carried out.
 */
typedef double (*msm4g_smoothing_handler)(double rho,int derivative);

/** @brief The collection of simulation parameters.
 *
 * The simulation parameters such as the cut-off threshold, the type 
 * of interpolating polynomials, the number of multi-levels are defined
 * in this structures.
 */
typedef struct SimulationParameters
{
    double                  a;                        /**< Cut-off parameter */
    msm4g_smoothing_handler msm4g_smoothing_function; /**< Smoothing function */
} SimulationParameters;

/** @brief The geomatric boundaries of the simulation.
 *
 * The algorithm requires that the simulation is confined into a rectangular 
 * box. The location and the width vectors are sufficient to desribe the domain.
 *
 * @todo If the width along on one of the axis is zero or too small, then
 * the bin should be allowed to extend beyond the simulation box.
 */
typedef struct SimulationBox
{
    D3Vector location; /**< Location of the lower left corner */
    D3Vector width;    /**< Width of each side */
} SimulationBox;

/** @brief Top level information about the simulation.
 *
 * This structure holds everything needed for the simulation takes part.
 */
typedef struct Simulation
{
    struct SimulationParameters simulationParameters; /**< Parameters of the algorithm. */
    struct SimulationBox        simulationBox;        /**< Geometry of the simulation. */
    struct LinkedList           *particles;           /**< The collection of the particles in the SimulationBox */
} Simulation;

#endif /* MSM4G_TYPES_H */
