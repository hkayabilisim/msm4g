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
    D3Vector fshort;        /**< force produced by short-range interactions. */
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
    struct SimulationParameters simulationParameters; /**< @brief All of the parameters of the algorithm. */
    struct SimulationBox        simulationBox;        /**< Geometry of the simulation. */
    struct LinkedList           *particles;           /**< The collection of the particles in the SimulationBox */
} Simulation;

/** @brief An abstract 3-dimensional grid structure.
 *
 * This AbstractGrid structure defines the interfaces that the grid implementations
 * should follow. Its main aim is to hide the details of different grid implementations
 * as much as possible and increase the code reuse.
 
 * Independant from the implementation, a 3D simulation domain is repeatedly divided into
 * grids whose nodes are stored is internal data structures of the implemention grid object
 * such as DenseGrid or SparseGrid. The visualiazation of an AbstractGrid is
 * shown below. In this 2D example (3D case is the straightforward generalization)
 * there are nx  and ny grid nodes in each axis and the lattice spacing is h.
 *
 @verbatim
 (ny-1)*h +=======+=======+=======+=======+
          |       |       |       |       |
 .        |       |       |       |       |
 .        |       |       |       |       |
 .        +=======+=======+=======+=======+
          |       |       |       |       |
          |       |       |       |       |
          |       |       |       |       |
 h        +=======+=======+=======+=======+
          |       |       |       |       |
          |       |       |       |       |
          |       |       |       |       |
 0        +=======+=======+=======+=======+
          0       h       2h     ....  (nx-1)*h
 
 @endverbatim
 *
 * In order to use AbstractGrid in the code, one should first use a constructor
 * designed for an implementation. After that, the methods msm4g_grid_xxx can be 
 * be used without dealing with details regarding the grid implementation.
 *
 * @code
   int nx=10; 
   int ny=10;
   int nz=10;
   double h=1.0;
   double element;
   AsbtractGrid *grid;
 
   grid = msm4g_dense_grid_new(nx,ny,nz,h);
   
   grid->reset(grid, 0.0);
 
   grid->setElement(grid,3,4,5,10.0);
 
   element = grid->getElement(grid,3,4,5);
 
   msm4g_grid_destroy(&grid);
 * @endcode
 */
typedef struct AbstractGrid
{
    /** @brief The lattice spacing.
     *
     * The distance between the grid nodes along each axis.
     * Please note that, one can also use different spacings along
     * each direction but simplicity, it is chosen the same for all directions.
     */
    double h;
    
    /** @brief The number of nodes along x-axis.
     *
     * This member defines the number of nodes
     * in the x-axis. Hence the number of boxes aling the axis will
     * nx-1. The coordinate of the last node will be (nx-1)*h.
     */
    int nx;
    
    /** @brief The number of nodes along y-axis.
     *
     * This member defines the number of nodes
     * in the y-axis. Hence the number of boxes aling the axis will
     * ny-1. The coordinate of the last node will be (ny-1)*h.
     */
    int ny;
    
    /** @brief The number of nodes along z-axis.
     *
     * This member defines the number of nodes
     * in the z-axis. Hence the number of boxes aling the axis will
     * nz-1. The coordinate of the last node will be (nz-1)*h.
     */
    int nz;
    
    /** @brief A pointer to the grid constructor.
     *
     * This function pointer determines the prototype of a grid
     * constructor which every grid implementations should create
     * a one.
     */
    struct AbstractGrid * (*constructor) (int nx,int ny,int nz, double h);
    
    /** @brief A pointer to the grid destructor.
     *
     * Every grid implementations should define a function satisfying 
     * this destructor prototype. But the implementations should only
     * deallocate their internal data structures.
     */
    void                  (*destructor)  (struct AbstractGrid **grid);


    /** @brief A pointer to the setter function.
     *
     * Every grid implementation should provide a setter function to set
     * a specific value to the given coordinates.
     */
    void   (*setElement)  (struct AbstractGrid *grid,int i,int j,int k,double value);
    
    /** @brief A pointer to the function prototype for getting a single element.
     *
     * The grid implementations should implement a function complying to this
     * prototype.
     */
    double (*getElement)   (struct AbstractGrid *grid,int i,int j,int k);
    
    /** @brief A pointer to the reset function.
     * 
     * Sets all elements in the grid to the given value.
     */
    void   (*reset)        (struct AbstractGrid *grid,double value);


} AbstractGrid;

/** @brief A DenseGrid implementation for the AbstractGrid.
 *
 * In this implementation, the nodes in the grid are stored 
 * consecutively in a dense matrix which is flattened into a
 * long vector with size nx * ny * nz.
 *
 * The constructor, destructor, setter and other methods 
 * related with the implementation are defined in msm4g_grid_dense_xxx
 * functions. They are later attached to the object during the 
 * construction of the object.
 */
typedef struct DenseGrid
{
    /** @brief The members inherited from AbstractGrid.
     *
     * The members common to all grid implementations are already
     * defined in the AbstractGrid structure. In order to get rid of
     * copy-paste of all members from AbstractGrid to DenseGrid, 
     * the content of AbstractGrid is inserted at the top of the 
     * DenseGrid structure. 
     *
     * @note It is very important to inlide the AbstractGrid at the
     * top because it is the only way to type cast a (DenseGrid *) to
     * (AbstractGrid *) or visa versa in a secure way.
     *
     * The members inherited from the AbstractGrid should not
     * meant to be directly accessed. Therefore the name will never be used
     * in the code. That's why, an obscure "_" symbol is used.
     */
    struct AbstractGrid _;
    
    /** @brief Internal data structure to store all the values on the nodes.
     *
     * This is the internal data structure to hold all of the
     * elements of the grid. It is important to define this member
     * after inherided members. The data array should be allocated in
     * the constructor of the grid and properly deallocated in the
     * destructor.
     */
    double *data;
    
} DenseGrid;

#endif /* MSM4G_TYPES_H */
