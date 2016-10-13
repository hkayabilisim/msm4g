#include "msm4g_lib.h"

int msm4g_linkedlist_size(LinkedList *list)
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

LinkedList *msm4g_linkedlist_new()
{
    LinkedList *newlist ;
    newlist = malloc(sizeof(LinkedList));
    newlist->head = NULL;
    newlist->tail = NULL;
    return newlist;
}

void msm4g_linkedlist_add(LinkedList *list,void *data)
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

void msm4g_linkedlist_destroy(LinkedList *list)
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

void msm4g_forwardeuler(double *r,double *v,double *a,double dt,int n,int d,double *m,double G)
{
    int i;
    msm4g_acceleration(a,r,n,d,m,G);
    for (i=0;i<n*d;i++)
        r[i] += v[i]*dt ;
    for (i=0;i<n*d;i++)
        v[i] += a[i]*dt ;
}

void msm4g_addbinlist(Bin **binlist,Body *body)
{
    
    
}

void msm4g_box_update(SimulationBox *box,LinkedList *list,double margin)
{
    const int DIM = 3;
    double min[DIM];
    double max[DIM];
    int i;
    LinkedListElement *curr ;
    Body *body;
    for (i=0;i<DIM;i++)
    {
        min[i] = DBL_MAX;
        max[i] = DBL_MIN;
    }
    
    curr = list->head;
    while (curr != NULL)
    {
        body = (Body *) (curr->data);
        for (i=0; i<DIM; i++)
        {
            if (body->r[i] >= max[i])
                max[i] = body->r[i];
            if (body->r[i] <= min[i])
                min[i] = body->r[i];
        }
        curr = curr->next;
    }
}

Body *msm4g_body_reset(Body *body)
{
    int i;
    body->m = 0.0;
    for (i=0; i<3; i++)
    {
        body->r[i] = 0.0;
        body->v[i] = 0.0;
        body->f[i] = 0.0;
    }
    return body;
}

Body *msm4g_body_rand(int n)
{
    Body *body;
    int i,j;
    body = malloc(sizeof(Body)*n);
    for (j=0;j<n;j++)
    {
        (body+j)->m = (double)rand()/RAND_MAX;
        for (i=0; i<3; i++)
        {
            (body+j)->r[i] = (double)rand()/RAND_MAX;
            (body+j)->v[i] = (double)rand()/RAND_MAX;
            (body+j)->f[i] = (double)rand()/RAND_MAX;
        }
    }
    return body;
}

Body *msm4g_body_new(double mass,double *location,double *velocity)
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

void msm4g_body_print(Body *body)
{
    printf("m:%8.2E ",body->m);
    printf("r:%8.2E %8.2E %8.2E ",body->r[0],body->r[1],body->r[2]);
    printf("v:%8.2E %8.2E %8.2E ",body->v[0],body->v[1],body->v[2]);
    printf("f:%8.2E %8.2E %8.2E ",body->f[0],body->f[1],body->f[2]);
    printf("\n");
    
}

void msm4g_acceleration(double *a,double *r,int n,int d,double *m,double G)
{
    int i,j,k;
    double rij2,rij3,acc;
    msm4g_zeros(a,n*d);
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


void msm4g_energy(double *pot,double *kin,double *tot,double *r,double *v,int n,int d,double *m,double G)
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



void msm4g_print_energy(double pot0,double kin0,double tot0,double pot,double kin,double tot)
{
    printf("[pot: %10.3E] [kin: %10.3E] [tot: %10.3E] [err: %10.3E]\n",pot,kin,tot,(tot-tot0)/tot0);
}

void msm4g_leapfrog(double *r,double *v,double *a,double *a1,double dt,int n,int d,double *m,double G)
{
    int i;
    msm4g_acceleration(a,r,n,d,m,G);
    for (i=0;i<n*d;i++)
        r[i] += v[i]*dt + a[i]*0.5*dt*dt;
    msm4g_acceleration(a1,r,n,d,m,G);
    for (i=0;i<n*d;i++)
        v[i] += (a[i]+a1[i])*0.5*dt;
}

void msm4g_zeros(double *x,int n)
{
    int i;
    for (i=0;i<n;i++)
        x[i]=0.0;
}
