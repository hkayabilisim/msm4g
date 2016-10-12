#include "nbody_lib.h"

/** @brief Running all unit tests sequentially. */
void nbody_unit_test_all()
{
    bool status;
    bool (*unitTests[4])() = {
        nbody_unit_test_1,
        nbody_unit_test_2,
        nbody_unit_test_3,
        nbody_unit_test_4
    };
    int testid, numberoftests,failed;
    
    numberoftests = sizeof(unitTests)/sizeof(unitTests[0]);
    failed=0;
    for (testid = 0; testid < numberoftests ; testid++)
    {
        status = unitTests[testid]();
        if (status == false)
        {
            failed++;
            fprintf(stderr,"Unit test %3d: %s\n",testid+1,"failed");
        }
    }
    printf("%d of %d passed, %d failed\n",numberoftests-failed,numberoftests,failed);
}

/** @brief Checking size() function of the linked list implementation. */
bool nbody_unit_test_1()
{
    LinkedList *list;
    Body x[10];
    int i;
    
    list = nbody_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        nbody_resetbody(x+i);
        nbody_linkedlist_add(list,&x[i]);
    }
    if (nbody_linkedlist_size(list) != 10)
    {
        nbody_linkedlist_destroy(list);
        return false;
    }
    nbody_linkedlist_destroy(list);
    return true;
}

/** @brief A test function for linked list implementation.
 *
 *
 */
bool nbody_unit_test_2()
{
    LinkedList *list;
    LinkedListElement *curr;
    Body x[10];
    Body *y;
    bool status = true;
    int i;
    
    list = nbody_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        x[i].m = i;
        nbody_linkedlist_add(list,&x[i]);
    }
    
    /**< Iterating from tail to head */
    i=9;
    curr=list->tail;
    while (curr != NULL)
    {
        y =  (Body *)curr->data;
        if (fabs(i - (y->m)) > DBL_EPSILON)
        {
            status = false;
            break;
        }
        curr = curr->prev;
        i--;
    }
    nbody_linkedlist_destroy(list);

    return status;
}

/** @brief Traversing from head to tail. */
bool nbody_unit_test_3()
{
    LinkedList *list;
    LinkedListElement *curr;
    Body x[10];
    Body *y;
    bool status = true;
    int i;
    
    list = nbody_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        x[i].m = i;
        nbody_linkedlist_add(list,&x[i]);
    }
    
    /**< Iterating from tail to head */
    i=0;
    curr=list->head;
    while (curr != NULL)
    {
        y =  (Body *)curr->data;
        if (fabs(i - (y->m)) > DBL_EPSILON)
        {
            status = false;
            break;
        }
        curr = curr->next;
        i++;
    }
    nbody_linkedlist_destroy(list);
    
    return status;
}

/** @brief Reading body configurations from a text file.
 * 
 * The masses, locations and velocities of three bodies are read from
 * eight.ini file and compared to the true values. If there is a slight
 * variation, the test returns with fail.
 *
 * @return true if the test is succesfull, false otherwise.
 */
bool nbody_unit_test_4()
{
    FILE *fp;
    double mass,r[3],v[3];
    Body bodies[3], bodiesInFile[3];
    int i,j;
    int ibody;
    
    for (i=0;i<3;i++)
    {
        nbody_resetbody(&bodies[i]);
        nbody_resetbody(&bodiesInFile[i]);
        bodies[i].m = 1.0;
    }
    bodies[0].r[0] =   0.97000436 ;
    bodies[0].r[1] =  -0.24308753 ;
    bodies[0].v[0] =   0.466203685 ;
    bodies[0].v[1] =   0.43236573 ;

    bodies[1].r[0] =  -0.97000436;
    bodies[1].r[1] =   0.24308753;
    bodies[1].v[0] =   0.466203685;
    bodies[1].v[1] =   0.43236573;
    
    bodies[2].v[0] =  -0.93240737;
    bodies[2].v[1] =  -0.86473146;

    fp = fopen("data/eight.ini","r");
    if (fp == NULL) return false;
    
    ibody = 0;
    while (true)
    {
        int ret = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&mass,&r[0],&r[1],&r[2],&v[0],&v[1],&v[2]);
        if (ret == 7)
        {
            bodiesInFile[ibody].m = mass;
            for (i=0;i<3;i++)
            {
                bodiesInFile[ibody].r[i] = r[i];
                bodiesInFile[ibody].v[i] = v[i];
            }
            ibody++;
        } else if (ret == EOF)
            break;
        else
        {
            fclose(fp);
            return false;
        }
            
    }
    
    for (i=0;i<3;i++) /* for each body */
    {
        if (fabs(bodies[i].m-bodiesInFile[i].m)>DBL_EPSILON )
            return false;
        for (j=0;j<3;j++) /* for each dimension */
        {
            if (fabs(bodies[i].r[j]-bodiesInFile[i].r[j]) > DBL_EPSILON)
                return false;
            if (fabs(bodies[i].v[j]-bodiesInFile[i].v[j]) > DBL_EPSILON)
                return false;
        }
    }
    
    fclose(fp);
    return true;
}
int nbody_linkedlist_size(LinkedList *list)
{
    int size = 0;
    LinkedListElement *curr;
    curr = list->head;
    while (curr != NULL)
    {
        size = size + 1;
        curr=curr->next;
    }
    return size;
}
LinkedList *nbody_linkedlist_new()
{
    LinkedList *newlist ;
    newlist = malloc(sizeof(LinkedList));
    newlist->head = NULL;
    newlist->tail = NULL;
    return newlist;
}
void nbody_linkedlist_add(LinkedList *list,void *data)
{
    LinkedListElement *item;
    
    item = malloc(sizeof(LinkedListElement));
    item->data = data;
    item->next = NULL;
    item->prev = NULL;
    
    /* the list is empty */
    if ( list->tail == NULL)
    {
        list->head = item;
        list->tail = item;
    } else
    {
        item->prev = list->tail;
        list->tail->next = item;
        list->tail = item;
    }
}

void nbody_linkedlist_destroy(LinkedList *list)
{
    LinkedListElement *current, *next;
    current = list->head;
    while (current != NULL)
    {
        next=current->next;
        free(current);
        current=next;
    }
    free(list);
}



void nbody_forwardeuler(double *r,double *v,double *a,double dt,int n,int d,double *m,double G)
{
    int i;
    nbody_acceleration(a,r,n,d,m,G);
    for (i=0;i<n*d;i++)
        r[i] += v[i]*dt ;
    for (i=0;i<n*d;i++)
        v[i] += a[i]*dt ;
}

void nbody_addbinlist(Bin **binlist,Body *body)
{
    
    
}




void nbody_createBox(SimulationBox *box,LinkedList *list)
{
    LinkedList *current;
    int i;
    double minmax[3][2];
    for (i=0;i<3;i++)
    {
        minmax[i][0] = DBL_MIN;
        minmax[i][1] = DBL_MAX;
    }
    current = list;
    while (current != NULL)
    {
        
    }
}

Body *nbody_resetbody(Body *body)
{
    int i;
    body->m = 1.0;
    for (i=0; i<3; i++)
    {
        body->r[i] = 0.0;
        body->v[i] = 0.0;
        body->f[i] = 0.0;
    }
    return body;
}

Body *nbody_newrandbody()
{
    Body *body;
    int i;
    body = malloc(sizeof(Body));
    body->m = (double)rand()/RAND_MAX;
    for (i=0; i<3; i++)
    {
        body->r[i] = (double)rand()/RAND_MAX;
        body->v[i] = (double)rand()/RAND_MAX;
        body->f[i] = (double)rand()/RAND_MAX;
    }
    return body;
}

Body *nbody_newbody(double mass,double *location,double *velocity)
{
    Body *body;
    body = malloc(sizeof(Body));
    body->m = mass;
    body->r[0] = location[0];
    body->r[1] = location[1];
    body->r[2] = location[2];
    body->v[0] = velocity[0];
    body->v[1] = velocity[1];
    body->v[2] = velocity[2];
    body->f[0] = 0.0;
    body->f[1] = 0.0;
    body->f[2] = 0.0;
    return body;
}

Body *nbody_newzerobody()
{
    Body *body;
    int i;
    body = malloc(sizeof(Body));
    body->m = 0.0;
    for (i = 0 ; i < 3; i++)
    {
        body->v[i] = i;
        body->r[i] = i;
        body->f[i] = i;
    }
    return body;
}

void nbody_printbody(Body *body)
{
    printf("m:%8.2E ",body->m);
    printf("r:%8.2E %8.2E %8.2E ",body->r[0],body->r[1],body->r[2]);
    printf("v:%8.2E %8.2E %8.2E ",body->v[0],body->v[1],body->v[2]);
    printf("f:%8.2E %8.2E %8.2E ",body->f[0],body->f[1],body->f[2]);
    printf("\n");
    
}

void nbody_acceleration(double *a,double *r,int n,int d,double *m,double G)
{
    int i,j,k;
    double rij2,rij3,acc;
    nbody_zeros(a,n*d);
    for (i=0;i<n-1;i++)
    {
        for (j=i+1;j<n;j++)
        {
            rij2=0.0;
            for (k=0;k<d;k++)
                rij2 +=(r[d*i+k]-r[d*j+k])*(r[d*i+k]-r[d*j+k]);
            rij3=sqrt(rij2)*rij2;
            for (k=0;k<d;k++)
            {
                acc=G*(r[d*i+k]-r[d*j+k])/rij3;
                a[d*j+k] += m[i]*acc;
                a[d*i+k] -= m[j]*acc;
            }
        }
    }
}


void nbody_energy(double *pot,double *kin,double *tot,double *r,double *v,int n,int d,double *m,double G)
{
    int i,j,k;
    double rij2,v2;
    *pot=0; *kin=0; *tot=0;
    for (i=0;i<n-1;i++)
    {
        for (j=i+1;j<n;j++)
        {
            rij2=0;
            for (k=0;k<d;k++)
            {
                rij2 += pow(r[d*i+k]-r[d*j+k],2);
            }
            *pot -= m[i]*m[j]/sqrt(rij2);
        }
    }
    *pot *= G;
    
    for (i=0;i<n;i++)
    {
        v2 = 0;
        for (k=0;k<d;k++)
        {
            v2 += v[d*i+k]*v[d*i+k];
        }
        *kin += m[i]*v2;
    }
    *kin *= 0.5;
    
    *tot = (*pot) + (*kin);
    
}

void nbody_print_body(double *r,double *v,double *a,int i)
{
    printf("%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
           r[i],r[i+1],r[i+2],v[i],v[i+1],v[i+2],a[i],a[i+1],a[i+2]);
}

void nbody_print_energy(double pot0,double kin0,double tot0,double pot,double kin,double tot)
{
    printf("[pot: %10.3E] [kin: %10.3E] [tot: %10.3E] [err: %10.3E]\n",pot,kin,tot,(tot-tot0)/tot0);
}

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
void nbody_leapfrog(double *r,double *v,double *a,double *a1,double dt,int n,int d,double *m,double G)
{
    int i;
    nbody_acceleration(a,r,n,d,m,G);
    for (i=0;i<n*d;i++)
        r[i] += v[i]*dt + a[i]*0.5*dt*dt;
    nbody_acceleration(a1,r,n,d,m,G);
    for (i=0;i<n*d;i++)
        v[i] += (a[i]+a1[i])*0.5*dt;
}

int nbody_test_eight()
{
    int n;                 /* number of particles */
    int d;                 /* dimension */
    double *r,*v,*a,*m;    /* location, velocity, acceleration and masses */
    double *a1;            /* extra acceleration needed in leapfrog */
    double t;              /* current time */
    double dt=0.01;     /* time step */
    double T=100;          /* duration of the integration */
    double G=1;            /* constant */
    double pot0,kin0,tot0; /* potential, kinetic, total energy */
    double pot,kin,tot;    /* potential, kinetic, total energy */
    struct timeval t_start,t_finish,t_elapsed;
    nbody_read_ini("data/eight.ini",&n,&d,&m,&r,&v);
    
    a =calloc(n*d,sizeof(double));
    a1=calloc(n*d,sizeof(double));
    
    nbody_energy(&pot0,&kin0,&tot0,r,v,n,d,m,G);
    nbody_energy(&pot, &kin, &tot, r,v,n,d,m,G);
    nbody_print_energy(pot0,kin0,tot0,pot,kin,tot);
    
    nbody_read_ini("data/eight.ini",&n,&d,&m,&r,&v);
    t=0.0;
    gettimeofday(&t_start,NULL);
    while (t<T)
    {
        /* nbody_forwardeuler(r,v,a,dt,n,d,m,G); */
        nbody_leapfrog(r,v,a,a1,dt,n,d,m,G);
        t += dt;
    }
    gettimeofday(&t_finish,NULL);
    timersub(&t_finish, &t_start, &t_elapsed);
    
    printf("%ld.%06ld ", (long int)t_elapsed.tv_sec, (long int)t_elapsed.tv_usec);
    
    nbody_energy(&pot,&kin,&tot,r,v,n,d,m,G);
    nbody_print_energy(pot0,kin0,tot0,pot,kin,tot);
    
    free(r);
    free(v);
    free(a);
    free(a1);
    free(m);
    return 0;
}

void nbody_zeros(double *x,int n)
{
    int i;
    for (i=0;i<n;i++)
        x[i]=0.0;
}
