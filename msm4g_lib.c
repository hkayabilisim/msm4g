/** @file msm4g_lib.c
 * @brief The definitions of all core functions of MSM4G package.
 */
#include "msm4g_lib.h"

double msm4g_smoothing_C1(double rho,int derivative)
{
    double rho2 ;
    rho2 = rho*rho;
    
    if (derivative == 1)
    {
        if (rho >= 1.0)
            return -1.0/rho2;
        else
            return -rho;
    } else
    {
        if (rho >= 1.0)
            return 1.0/rho;
        else
            return 0.5*(3.0 - rho2);
    }
}

double msm4g_smoothing_C2(double rho,int derivative)
{
    double rho2,rho3,rho4;
    rho2 = rho*rho;

    if (derivative == 1)
    {
        if (rho >= 1.0)
            return -1.0/rho2;
        else
        {
            rho3 = rho2*rho;
            return 0.5*(-5.0*rho + 3.0*rho3);
        }

    } else
    {
    if (rho >= 1.0)
        return 1.0/rho;
    else
    {
        rho2 = rho*rho;
        rho4 = rho2*rho2;
        return 0.125*(15.0 - 10.0*rho2 + 3.0*rho4);
    }
    }
}

double msm4g_smoothing_C3(double rho,int derivative)
{
    double rho2,rho3,rho4,rho5,rho6;
    rho2 = rho*rho;
    if (derivative == 1)
    {
        if (rho >= 1.0)
            return -1.0/rho2;
        else
        {
            rho3 = rho2*rho;
            rho5 = rho2*rho3;
            return 0.125*(-35.0*rho + 42.0*rho3 - 15.0*rho5);
        }
    }
    else
    {
        if (rho >= 1.0)
            return 1.0/rho;
        else
        {
            rho2 = rho*rho;
            rho4 = rho2*rho2;
            rho6 = rho2*rho4;
            return 0.0625*(35.0 - 35.0*rho2 + 21*rho4 - 5*rho6);
        }
    }
}

void msm4g_force_short(LinkedList *binlist,double threshold, msm4g_smoothing_handler smoothing_function)
{
    Bin *bin;
    Bin *neighbor;
    LinkedListElement *currBin;
    LinkedListElement *neighborBin;
    
    currBin=binlist->head;
    while (currBin != NULL)
    {
        bin = (Bin *)currBin->data;
        msm4g_force_short_withinBin(bin->particles,threshold,smoothing_function);
        
        neighborBin = bin->neighbors->head;
        while (neighborBin != NULL)
        {
            neighbor = (Bin *)neighborBin->data;
            if (neighbor->cantorindex > bin->cantorindex )
            {
                printf("Source bin index: "); msm4g_i3vector_print(&(bin->index)); printf("\n");
                printf("target bin index: "); msm4g_i3vector_print(&(neighbor->index)); printf("\n");
            }
            neighborBin=neighborBin->next;
        }
        currBin = currBin->next;
    }
}

void msm4g_force_short_withinBin(LinkedList *particles, double threshold, msm4g_smoothing_handler smoothing_function)
{
    Particle *particleI;
    Particle *particleJ;
    LinkedListElement *currI;
    LinkedListElement *currJ;
    

    
    
    currI = particles->head;
    while (currI->next != NULL)
    {
        particleI = (Particle *)currI->data;
        currJ = currI->next;
        while (currJ != NULL)
        {
            particleJ = (Particle *)currJ->data;
            printf("%d versus %d\n",particleI->index,particleJ->index);
            msm4g_force_short_particlePair(particleI,particleJ,threshold,smoothing_function);
            
            currJ = currJ->next;
        }
        currI = currI->next;
    }
}

void msm4g_force_short_particlePair(Particle *particleI, Particle *particleJ, double a, msm4g_smoothing_handler gamma)
{
    D3Vector rij;
    double r,r2;
    double a2;
    double smoothedkernel, smoothedkernelderivative;
    
    /* rij = rj - ri */
    msm4g_d3vector_daxpy(&rij, &(particleJ->r), -1.0, &(particleI->r));
    /* r^2 = |rij|^2 */
    r2=msm4g_d3vector_normsquare(&rij);
    r=sqrt(r2);
    a2=a*a;
    /* If only the particles are closer than cut-off distance */
    if (r2 < a2)
    {
        smoothedkernel = 1.0/r - gamma(r/a,0)/a;
        smoothedkernelderivative = 1.0/r2 - gamma(r/a,1)/a2;
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
    return sqrt(msm4g_d3vector_normsquare(x));
}

double msm4g_d3vector_normsquare(D3Vector *x)
{
    double normsquare = 0.0;
    int i;
    for (i=0; i<3; i++)
    {
        normsquare += x->value[i] * x->value[i];
    }
    return normsquare;
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
    Particle *particle;
    for (i=0;i<3;i++)
    {
        min[i] = DBL_MAX;
        max[i] = DBL_MIN;
    }
    
    curr = list->head;
    while (curr != NULL)
    {
        particle = (Particle *) (curr->data);
        for (i=0; i<3; i++)
        {
            if (particle->r.value[i] >= max[i])
                max[i] = particle->r.value[i];
            if (particle->r.value[i] <= min[i])
                min[i] = particle->r.value[i];
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

void msm4g_box_translate(SimulationBox *box, LinkedList *particles,D3Vector delta)
{
    int i;
    int n;
    n=3;
    LinkedListElement *curr;
    Particle *particle;
    
    for (i=0;i<n;i++)
        box->location.value[i] += delta.value[i];
    curr = particles->head;
    while (curr != NULL)
    {
        particle = (Particle *)curr->data;
        for (i=0;i<n;i++)
            particle->r.value[i] += delta.value[i];
        curr = curr->next;
    }
}

void msm4g_box_translateToOrigin(SimulationBox *box, LinkedList *particles)
{
    D3Vector delta;
    delta.value[0] =  -box->location.value[0];
    delta.value[1] =  -box->location.value[1];
    delta.value[2] =  -box->location.value[2];
    msm4g_box_translate(box,particles,delta);
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
    bin->particles = msm4g_linkedlist_new();
    bin->neighbors = msm4g_linkedlist_new();
    msm4g_i3vector_set(&(bin->index), index.value[0], index.value[1], index.value[2]);
    bin->cantorindex = msm4g_math_cantor(index.value, 3);
    return bin;
}

LinkedList *msm4g_bin_generate(SimulationBox *box,LinkedList *particles,double binwidth)
{
    LinkedList *binlist;
    LinkedListElement *curr;
    Particle *particle;
    Bin *bin;
    I3Vector binindex;
    int i;

    binlist = msm4g_linkedlist_new();
    
    curr = particles->head;
    while (curr != NULL)
    {
        particle = (Particle *)curr->data;
        msm4g_i3vector_set(&binindex, -1, -1, -1);
        for (i=0;i<3;i++)
        {
            binindex.value[i] = floor(particle->r.value[i]/binwidth);
        }
        bin = msm4g_bin_searchByIndex(binlist,binindex);
        if (bin == NULL)
        {
            bin = msm4g_bin_new(binindex);
            msm4g_linkedlist_add(binlist, bin);
        }
        msm4g_linkedlist_add(bin->particles, particle);
        curr = curr->next;
    }
    
    msm4g_bin_findneighbors(binlist);
    
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
    Particle *particle;
    LinkedListElement *curr;
    
    printf("[%d,%d,%d] cantor index: %d\n",
           bin->index.value[0],bin->index.value[1],bin->index.value[2],
           bin->cantorindex);
    
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
    
    /* Printing particles */
    curr = bin->particles->head;
    while (curr != NULL)
    {
        particle = (Particle *)curr->data;
        printf("  particle: %d [%f,%f,%f]\n",particle->index,
               particle->r.value[0],
               particle->r.value[1],
               particle->r.value[2]);
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
        msm4g_linkedlist_destroy(bin->particles);
        msm4g_linkedlist_destroy(bin->neighbors);
        curr = curr->next;
    }
    msm4g_linkedlist_destroyWithData(binlist);
}

Particle *msm4g_particle_reset(Particle *particle)
{
    particle->m = 0.0;
    msm4g_d3vector_set(&(particle->r), 0, 0, 0);
    msm4g_d3vector_set(&(particle->v), 0, 0, 0);
    msm4g_d3vector_set(&(particle->f), 0, 0, 0);
    return particle;
}

Particle **msm4g_particle_rand(int n)
{
    Particle **particle;
    int i,j;
    
    particle = malloc(sizeof(Particle *)*n);
    
    for (j=0;j<n;j++)
    {
        particle[j] = msm4g_particle_empty();
        particle[j]->m = (double)rand()/RAND_MAX;
        for (i=0; i<3; i++)
        {
            particle[j]->r.value[i] = (double)rand()/RAND_MAX;
            particle[j]->v.value[i] = (double)rand()/RAND_MAX;
            particle[j]->f.value[i] = (double)rand()/RAND_MAX;
        }
    }
    return particle;
}

LinkedList *msm4g_particle_read(const char *filename)
{
    LinkedList *particles;
    Particle *particle;
    FILE *fp;
    int iparticle;
    double mass;
    double r[3];
    double v[3];
    
    particles = msm4g_linkedlist_new();
    fp = fopen(filename,"r");
    if (fp == NULL)
    {
        return particles;
    }
    
    iparticle = 0;
    while (true)
    {
        int ret = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&mass,&r[0],&r[1],&r[2],&v[0],&v[1],&v[2]);
        if (ret == 7)
        {
            particle = msm4g_particle_new(mass, r, v);
            msm4g_linkedlist_add(particles, particle);
            iparticle++;
        } else if (ret == EOF)
            break;
    }
    fclose(fp);
    
    return particles;
}

Particle *msm4g_particle_empty()
{
    static int index = 0;
    int i;
    Particle *particle;
    
    particle = malloc(sizeof(Particle));
    particle->index = index;
    particle->m = 0.0;
    for (i=0; i<3; i++)
    {
        particle->r.value[i] = 0.0;
        particle->v.value[i] = 0.0;
        particle->f.value[i] = 0.0;
    }
    
    index++;
    return particle;
}

Particle *msm4g_particle_new(double mass,double *location,double *velocity)
{
    Particle *particle;
    particle = msm4g_particle_empty();
    particle->m = mass;
    particle->r.value[0] = location[0];
    particle->r.value[1] = location[1];
    particle->r.value[2] = location[2];
    particle->v.value[0] = velocity[0];
    particle->v.value[1] = velocity[1];
    particle->v.value[2] = velocity[2];
    return particle;
}

void msm4g_particle_print(Particle *particle)
{
    printf("i:%d ",particle->index);
    printf("m:%8.2E ",particle->m);
    printf("r:%8.2E %8.2E %8.2E ",particle->r.value[0],particle->r.value[1],particle->r.value[2]);
    printf("v:%8.2E %8.2E %8.2E ",particle->v.value[0],particle->v.value[1],particle->v.value[2]);
    printf("f:%8.2E %8.2E %8.2E ",particle->f.value[0],particle->f.value[1],particle->f.value[2]);
    printf("\n");
    
}

void msm4g_particle_printlist(LinkedList *particlelist)
{
    Particle *particle;
    LinkedListElement *curr;
    int i=0;
    curr = particlelist->head;
    while (curr != NULL)
    {
        particle = curr->data;
        printf("Particle index: %d\n",i);
        msm4g_particle_print(particle);
        
        i++;
        curr=curr->next;
    }
}

void msm4g_particle_destroy(Particle *particle)
{
    free(particle);
    particle=NULL;
}

void msm4g_particle_destroyarray(Particle **particlearray,int length)
{
    int i;
    for (i=0; i<length; i++)
    {
        free(particlearray[i]);
        particlearray[i] = NULL;
    }
    free(particlearray);
    particlearray = NULL;
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

int msm4g_math_cantor(int *x,int n)
{
    int i,j,h;
    int hsum = 0;
    int jprod = 1;
    int isum=0;
    for (h=1; h<=n; h++)
    {
        jprod = 1;
        for (j=0; j<= h-1; j++)
        {
            isum = 0;
            for (i=1; i<=h; i++)
            {
                isum += x[i-1];
            }
            jprod *= isum + j;
        }
        hsum += jprod/msm4g_math_factorial(h);
    }
    return hsum;
}

int msm4g_math_factorial(int n)
{
    int factorial = 1;
    int i;
    if (n==0) return 1;
    for (i=1; i<= n; i++)
    {
        factorial *= i;
    }
    return factorial;
}

