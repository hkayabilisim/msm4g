/** Library.
 *
 */

#include "nbody_lib.h"


void nbody_linkedlist_test1()
{
    LinkedList *list; /** The list */
    LinkedListElement *curr;
    Body x[10];
    int i;
    
    list = nbody_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        nbody_resetbody(x+i);
        nbody_linkedlist_add(list,&x[i]);
    }
    printf("Size of the list: %d\n",nbody_linkedlist_size(list));

    printf("Printfing from tail to head\n");
    curr=list->tail;
    while (curr != NULL)
    {
        nbody_printbody((Body *)curr->data);
        curr = curr->prev;
    }
    
    nbody_linkedlist_destroy(list);
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

void nbody_freebodylist(BodyList **bodylist)
{
    BodyList *current,*next;
    current = *bodylist ;
    while (current->next != NULL)
    {
        next = current->next;
        free(current->body);
        free(current);
        current = next;
    }
    free(current->body);
    free(current);
    *bodylist = NULL;
}

void nbody_addbinlist(Bin **binlist,Body *body)
{
    
    
}

BodyList *nbody_addbodylist(BodyList **bodylist,Body *body)
{
    BodyList *newnode ;
    BodyList *current;
    
    newnode = malloc(sizeof(BodyList));
    newnode->body = body;
    newnode->next = NULL;
    
    current = *bodylist;
    if (current == NULL)
    {
        *bodylist = newnode;
    } else
    {
        while (current->next != NULL)
        {
            current = current->next;
        }
        current->next = newnode;
    }
    
    return newnode;
}

void nbody_printbodylist(BodyList *head)
{
    BodyList *current;
    current = head;
    while (current != NULL)
    {
        nbody_printbody(current->body);
        current = current->next;
    }
}

typedef struct node * pNodeStruct;
typedef struct node {
    int data;
    pNodeStruct pNext;
} NodeStruct;

void nbody_createBox(SimulationBox *box,BodyList *list)
{
    BodyList *current;
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
    body->m = 0.5;
    for (i=0; i<3; i++)
    {
        body->r[i] = 1.0;
        body->v[i] = 2.0;
        body->f[i] = 3.0;
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

void nbody_read_ini(char *filename,int *n,int *d,double **m,double **r,double **v)
{
    int i,j;
    FILE *fp;
    fp = fopen(filename,"r");
    fscanf(fp,"%d %d",n,d);
    *r =calloc((*n)*(*d),sizeof(double));
    *m =calloc(*n,sizeof(double));
    *v =calloc((*n)*(*d),sizeof(double));
    for (i=0;i< (*n); i++) fscanf(fp,"%lf",*m+i);
    for (i=0;i< (*n);i++)
    {
        for (j=0;j< (*d);j++) fscanf(fp,"%lf",*r + (*d)*i+j);
        for (j=0;j< (*d);j++) fscanf(fp,"%lf",*v + (*d)*i+j);
    }
    fclose(fp);
    
}

int nbody_read_bodies(char *filename,int *n,BodyList **bodylist)
{
    int i,j;
    int d;
    double *m,*r,*v;
    FILE *fp;
    fp = fopen(filename,"r");
    fscanf(fp,"%d %d",n,&d);
    
    if (d != 3)
    {
        fprintf(stderr,"The simulation should be in 3D");
        return NBODY_ERR_INPUT_ISNOT_3D;
    }
    
    m =calloc(*n,sizeof(double));
    r =calloc((*n)*3,sizeof(double));
    v =calloc((*n)*3,sizeof(double));
    for (i=0;i< (*n); i++) fscanf(fp,"%lf",m+i);
    for (i=0;i< (*n);i++)
    {
        for (j=0;j< 3;j++) fscanf(fp,"%lf",r + (3)*i+j);
        for (j=0;j< 3;j++) fscanf(fp,"%lf",v + (3)*i+j);
    }
    fclose(fp);
    
    for (i=0; i < *n ;i++)
    {
        Body *body = nbody_newbody(m[i], r+3*i, v+3*i);
        //printf("Adding body\n");
        //nbody_printbody(body);
        
        nbody_addbodylist(bodylist, body);
    }
    
    free(m);
    free(r);
    free(v);
    return 0;
}

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
