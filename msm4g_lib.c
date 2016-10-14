/** @file msm4g_lib.c
 * @brief The definitions of all core functions of MSM4G package.
 */
#include "msm4g_lib.h"

void msm4g_force_short(LinkedList *binlist,double threshold)
{
    Bin *bin;
    Bin *binNeighbor;
    LinkedListElement *curr;
    LinkedListElement *currNeighbor;
    curr=binlist->head;
    while (curr != NULL)
    {
        bin = (Bin *)curr->data;
        
        currNeighbor = bin->neighbors->head;
        while (currNeighbor != NULL)
        {
            binNeighbor = (Bin *)currNeighbor->data;
            currNeighbor=currNeighbor->next;
            printf("Source bin index: "); msm4g_i3vector_print(&(bin->index)); printf("\n");
            printf("target bin index: "); msm4g_i3vector_print(&(binNeighbor->index)); printf("\n");
            
            

        }
        
        
        curr = curr->next;
    }
}
void msm4g_d3vector_set(D3Vector *d3vector,double x,double y,double z)
{
    d3vector->value[0] = x;
    d3vector->value[1] = y;
    d3vector->value[2] = z;
}

void msm4g_d3vector_daxpy(D3Vector *z,D3Vector *x,double a,D3Vector *y)
{
    int i;
    for (i=0; i<3; i++)
    {
        z->value[i] = x->value[i] + a*y->value[i];
    }
}

double msm4g_d3vector_norm(D3Vector *x)
{
    double norm = 0.0;
    int i;
    for (i=0; i<3; i++)
    {
        norm += x->value[i] * x->value[i];
    }
    norm = sqrt(norm);
    return norm;
}

void msm4g_d3vector_print(D3Vector *x)
{
    int i;
    for (i=0; i<3; i++)
    {
        printf("%10.3E ",x->value[i]);
    }
}

void msm4g_i3vector_set(I3Vector *i3vector,int x,int y,int z)
{
    i3vector->value[0] = x;
    i3vector->value[1] = y;
    i3vector->value[2] = z;
}

void msm4g_i3vector_copy(I3Vector *to, I3Vector from)
{
    to->value[0] = from.value[0];
    to->value[1] = from.value[1];
    to->value[2] = from.value[2];
}

Boolean msm4g_i3vector_isequal(I3Vector *x,I3Vector *y)
{
    int i;
    if (x == NULL || y == NULL)
        return false;
    for (i=0;i<3;i++)
    {
        if (x->value[i] != y->value[i])
            return false;
    }
    return true;
}

void msm4g_i3vector_print(I3Vector *x)
{
    printf("[%d,%d,%d] ",x->value[0],x->value[1],x->value[2]);
}

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

void *msm4g_linkedlist_get(LinkedList *list, int index)
{
    void *data = NULL;
    LinkedListElement *curr;
    int i;
    
    curr = list->head;
    i=0;
    while (curr != NULL)
    {
        if (i == index)
        {
            data = curr->data;
            break;
        }
        i++;
        curr = curr->next;        
    }
    return data;
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
    list = NULL;
}

void msm4g_linkedlist_destroyWithData(LinkedList *list)
{
    LinkedListElement *curr;
    curr = list->head;
    while (curr != NULL)
    {
        free(curr->data);
        curr = curr->next;
    }
    
    msm4g_linkedlist_destroy(list);
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

SimulationBox *msm4g_box_new()
{
    SimulationBox *box;
    box = malloc(sizeof(SimulationBox));
    msm4g_d3vector_set(&(box->location), 0.0,0.0,0.0);
    msm4g_d3vector_set(&(box->width), 1.0,1.0,1.0);
    return box;
}

void msm4g_box_update(SimulationBox *box,LinkedList *list,double margin)
{
    double min[3];
    double max[3];
    double oldwidth,newwidth;
    int i;
    LinkedListElement *curr ;
    Body *body;
    for (i=0;i<3;i++)
    {
        min[i] = DBL_MAX;
        max[i] = DBL_MIN;
    }
    
    curr = list->head;
    while (curr != NULL)
    {
        body = (Body *) (curr->data);
        for (i=0; i<3; i++)
        {
            if (body->r[i] >= max[i])
                max[i] = body->r[i];
            if (body->r[i] <= min[i])
                min[i] = body->r[i];
        }
        curr = curr->next;
    }
    
    for (i=0;i<3;i++)
    {
        box->location.value[i] = min[i];
        box->width.value[i] = max[i]-min[i];
    }
    
    if (margin > 0.0)
    {
        for (i=0;i<3;i++)
        {
            oldwidth = box->width.value[i];
            newwidth = oldwidth*(1.0+margin);
            box->location.value[i] -= 0.5*(newwidth-oldwidth);
            box->width.value[i] = newwidth;
        }
    }
}

void msm4g_box_translate(SimulationBox *box, LinkedList *bodies,D3Vector delta)
{
    int i;
    int n;
    n=3;
    LinkedListElement *curr;
    Body *body;
    
    for (i=0;i<n;i++)
        box->location.value[i] += delta.value[i];
    curr = bodies->head;
    while (curr != NULL)
    {
        body = (Body *)curr->data;
        for (i=0;i<n;i++)
            body->r[i] += delta.value[i];
        curr = curr->next;
    }
}

void msm4g_box_translateToOrigin(SimulationBox *box, LinkedList *bodies)
{
    D3Vector delta;
    delta.value[0] =  -box->location.value[0];
    delta.value[1] =  -box->location.value[1];
    delta.value[2] =  -box->location.value[2];
    msm4g_box_translate(box,bodies,delta);
}

void msm4g_box_print(SimulationBox *box)
{
    printf("Box Location: %f %f %f\n",
           box->location.value[0],
           box->location.value[1],
           box->location.value[2]);
    printf("Box Width   : %f %f %f\n",
           box->width.value[0],
           box->width.value[1],
           box->width.value[2]);
}

void msm4g_box_destroy(SimulationBox *box)
{
    free(box);
}

Bin *msm4g_bin_new(I3Vector index)
{
    Bin *bin;
    bin = malloc(sizeof(Bin));
    bin->bodies = msm4g_linkedlist_new();
    bin->neighbors = msm4g_linkedlist_new();
    msm4g_i3vector_set(&(bin->index), index.value[0], index.value[1], index.value[2]);
    return bin;
}

LinkedList *msm4g_bin_generate(SimulationBox *box,LinkedList *bodies,double binwidth)
{
    LinkedList *binlist;
    LinkedListElement *curr;
    Body *body;
    Bin *bin;
    I3Vector binindex;
    int i;

    binlist = msm4g_linkedlist_new();
    
    curr = bodies->head;
    while (curr != NULL)
    {
        body = (Body *)curr->data;
        msm4g_i3vector_set(&binindex, -1, -1, -1);
        for (i=0;i<3;i++)
        {
            binindex.value[i] = floor(body->r[i]/binwidth);
        }
        bin = msm4g_bin_searchByIndex(binlist,binindex);
        if (bin == NULL)
        {
            bin = msm4g_bin_new(binindex);
            msm4g_linkedlist_add(binlist, bin);
        }
        msm4g_linkedlist_add(bin->bodies, body);
        curr = curr->next;
    }
    
    msm4g_bin_findneighbors(binlist);
    
    curr = binlist->head;
    while (curr != NULL)
    {
        bin = (Bin *)curr->data;
        msm4g_bin_print(bin);
        curr = curr->next;
    }
    return binlist;
}

void msm4g_bin_findneighbors(LinkedList *binlist)
{
    Bin *bin;
    Bin *neighborBin;
    LinkedListElement *curr;
    I3Vector index;
    int i,j,k;
    
    curr = binlist->head;
    while (curr != NULL)
    {
        bin = (Bin *)curr->data;
        for (i=-1; i < 2; i++)
        {
            for (j=-1; j<2; j++)
            {
                for (k=-1; k<2 ; k++)
                {
                    if (i==0 && j==0 && k==0)
                        continue;
                    msm4g_i3vector_copy(&index, bin->index);
                    index.value[0] += i;
                    index.value[1] += j;
                    index.value[2] += k;
                    neighborBin = msm4g_bin_searchByIndex(binlist, index);
                    if (neighborBin != NULL)
                    {
                        msm4g_linkedlist_add(bin->neighbors,neighborBin);
                    }
                }
            }
        }
        curr = curr->next;
    }
}

Bin *msm4g_bin_searchByIndex(LinkedList *binlist,I3Vector index)
{
    LinkedListElement *curr;
    Bin *bin = NULL;
    curr = binlist->head;
    while (curr != NULL)
    {
        bin = (Bin *)curr->data;
        if (msm4g_i3vector_isequal(&index, &(bin->index)))
        {
           /*  printf("%d-%d %d-%d %d-%d\n",
                   index.value[0],bin->index.value[0],
                   index.value[1],bin->index.value[1],
                   index.value[2],bin->index.value[2]); */
            return bin;
        }
        curr = curr->next;
    }
    
    return NULL;
}

void msm4g_bin_print(Bin *bin)
{
    Bin *neighborBin;
    Body *body;
    LinkedListElement *curr;
    
    printf("[%d,%d,%d]\n",bin->index.value[0],bin->index.value[1],bin->index.value[2]);
    
    /* Printing neighbor bins */
    curr = bin->neighbors->head;
    while (curr != NULL)
    {
        neighborBin = (Bin *)curr->data;
        printf("  neighbor: [%d,%d,%d]\n",neighborBin->index.value[0],
                               neighborBin->index.value[1],
                               neighborBin->index.value[2]);

        curr = curr->next;
    }
    
    /* Printing bodies */
    curr = bin->bodies->head;
    while (curr != NULL)
    {
        body = (Body *)curr->data;
        printf("  body: %d [%f,%f,%f]\n",body->index,body->r[0],body->r[1],body->r[2]);
        curr = curr->next;
    }
}

void msm4g_bin_printlist(LinkedList *binlist)
{
    
}

void msm4g_bin_destroy(LinkedList *binlist)
{
    LinkedListElement *curr;
    Bin *bin;
    
    curr = binlist->head;
    while (curr != NULL)
    {
        bin = (Bin *)curr->data;
        msm4g_linkedlist_destroy(bin->bodies);
        msm4g_linkedlist_destroy(bin->neighbors);
        curr = curr->next;
    }
    msm4g_linkedlist_destroyWithData(binlist);
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

Body **msm4g_body_rand(int n)
{
    Body **body;
    int i,j;
    
    body = malloc(sizeof(Body *)*n);
    
    for (j=0;j<n;j++)
    {
        body[j] = msm4g_body_empty();
        body[j]->m = (double)rand()/RAND_MAX;
        for (i=0; i<3; i++)
        {
            body[j]->r[i] = (double)rand()/RAND_MAX;
            body[j]->v[i] = (double)rand()/RAND_MAX;
            body[j]->f[i] = (double)rand()/RAND_MAX;
        }
    }
    return body;
}

LinkedList *msm4g_body_read(const char *filename)
{
    LinkedList *bodies;
    Body *body;
    FILE *fp;
    int ibody;
    double mass;
    double r[3];
    double v[3];
    
    bodies = msm4g_linkedlist_new();
    fp = fopen(filename,"r");
    if (fp == NULL)
    {
        return bodies;
    }
    
    ibody = 0;
    while (true)
    {
        int ret = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&mass,&r[0],&r[1],&r[2],&v[0],&v[1],&v[2]);
        if (ret == 7)
        {
            body = msm4g_body_new(mass, r, v);
            msm4g_linkedlist_add(bodies, body);
            ibody++;
        } else if (ret == EOF)
            break;
    }
    fclose(fp);
    
    return bodies;
}

Body *msm4g_body_empty()
{
    static int index = 0;
    int i;
    Body *body;
    
    body = malloc(sizeof(Body));
    body->index = index;
    body->m = 0.0;
    for (i=0; i<3; i++)
    {
        body->r[i] = 0.0;
        body->v[i] = 0.0;
        body->f[i] = 0.0;
    }
    
    index++;
    return body;
}

Body *msm4g_body_new(double mass,double *location,double *velocity)
{
    Body *body;
    body = msm4g_body_empty();
    body->m = mass;
    body->r[0] = location[0];
    body->r[1] = location[1];
    body->r[2] = location[2];
    body->v[0] = velocity[0];
    body->v[1] = velocity[1];
    body->v[2] = velocity[2];
    return body;
}

void msm4g_body_print(Body *body)
{
    printf("i:%d ",body->index);
    printf("m:%8.2E ",body->m);
    printf("r:%8.2E %8.2E %8.2E ",body->r[0],body->r[1],body->r[2]);
    printf("v:%8.2E %8.2E %8.2E ",body->v[0],body->v[1],body->v[2]);
    printf("f:%8.2E %8.2E %8.2E ",body->f[0],body->f[1],body->f[2]);
    printf("\n");
    
}

void msm4g_body_printlist(LinkedList *bodylist)
{
    Body *body;
    LinkedListElement *curr;
    int i=0;
    curr = bodylist->head;
    while (curr != NULL)
    {
        body = curr->data;
        printf("Body index: %d\n",i);
        msm4g_body_print(body);
        
        i++;
        curr=curr->next;
    }
}

void msm4g_body_destroy(Body *body)
{
    free(body);
    body=NULL;
}

void msm4g_body_destroyarray(Body **bodyarray,int length)
{
    int i;
    for (i=0; i<length; i++)
    {
        free(bodyarray[i]);
        bodyarray[i] = NULL;
    }
    free(bodyarray);
    bodyarray = NULL;
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