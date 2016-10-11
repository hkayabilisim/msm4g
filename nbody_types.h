/** @brief Data structures
 *
 * The file contains all the
 */

/** TRUE corresponds to 1 */
#define TRUE  1
/** FALSE corresponds to 0 */
#define FALSE 0


typedef struct LinkedListElement
{
    void *data;
    struct LinkedListElement *next;
    struct LinkedListElement *prev;
} LinkedListElement;

typedef struct LinkedList
{
    struct LinkedListElement *head;
    struct LinkedListElement *tail;
} LinkedList;

typedef double Vector[3];

typedef struct Body
{
    double m;           /**< mass      */
    double r[3];        /**< location  */
    double v[3];        /**< velocity  */
    double f[3];        /**< force     */
} Body;

typedef struct BodyList
{
    struct Body *body;
    struct BodyList *next;
} BodyList;

typedef struct Bin
{
    int i;
    int j;
    int k;
    struct BodyList *bodylist;
    struct Bin *neighbors;
} Bin;

typedef struct SimulationParameters
{
    double a; /* Cut-off parameter */
} SimulationParameters;

typedef struct SimulationBox
{
    double location[3];
    double width[3];
} SimulationBox;

typedef struct Simulation
{
    struct SimulationParameters simulationParameters;
    struct SimulationBox        simulationBox;
} Simulation;


typedef struct Cell
{
    int i;
} Cell;
