/** @file msm4g_types.h
 * @brief The data type definitions required in the package.
 *
 * It is already included in @link msm4g_lib.h, so that it
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
 * represented with this simple Body structure containing basic features
 * such as mass, location, velocity, and force.
 *
 * @warning The index is supposed to be read-only in the user-end.
 * This variable is created when the body is allocated for the very first 
 * time inside the msm4g_body_empty() function. This means that user should
 * not allocate a body manually. He/She should always use msm4g_body functions
 * to handle the bodies.
 */
typedef struct Body
{
    int    index;       /**< Index of the body (see warning in the description) */
    double m;           /**< mass      */
    double r[3];        /**< location  */
    double v[3];        /**< velocity  */
    double f[3];        /**< force     */
} Body;

/** @brief The cubical domains obtained by dividing the simulation box.
 * 
 * Bin is a rectangular box created by dividing the simulation box.
 */
typedef struct Bin
{
    I3Vector index;         /**< The index of the bin balong each axis */
    LinkedList *neighbors;  /**< The list of neightbor bins containing at leeast one body */
    LinkedList *bodies;     /**< The list of bodies in the bin */
} Bin;

/** @brief The collection of simulation parameters.
 *
 * The simulation parameters such as the cut-off threshold, the type 
 * of interpolating polynomials, the number of multi-levels are defined
 * in this structures.
 */
typedef struct SimulationParameters
{
    double a; /**< Cut-off parameter */
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
    struct LinkedList          *bodies;               /**< The collection of the bodies in the SimulationBox */
} Simulation;

#endif /* MSM4G_TYPES_H */
