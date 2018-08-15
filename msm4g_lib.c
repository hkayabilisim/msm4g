/** @file msm4g_lib.c
 * @brief The definitions of all core functions of MSM4G package.
 */
#include "msm4g_lib.h"
#include "msm4g_bases.h"

void msm4g_anterpolation(AbstractGrid *gridmass,SimulationBox *box,LinkedList *particles,const BaseFunction *base)
{
    int i0,j0,k0;
    int i,j,k;
    double g,gold;
    double mass;
    int p = base->p;
    int v;
    double x,y,z;
    double x0,y0,z0;
    double rx_hx,ry_hy,rz_hz;
    int s_edge;
    double h = gridmass->h;
    double tx,ty,tz;
    double phix[MAX_POLY_DEGREE], phiy[MAX_POLY_DEGREE], phiz[MAX_POLY_DEGREE];
    LinkedListElement *curr;
    Particle *particle;
    
    x0 = box->location.value[0];
    y0 = box->location.value[1];
    z0 = box->location.value[2];

    s_edge = p/2-1;
    curr = particles->head;
    while (curr != NULL)
    {
        particle = (Particle *)curr->data;
        x =  particle->r.value[0];
        y =  particle->r.value[1];
        z =  particle->r.value[2];
        mass = particle->m;
        
        rx_hx = (x-x0) / h;
        ry_hy = (y-y0) / h;
        rz_hz = (z-z0) / h;
        
        i0 = floor(rx_hx) - s_edge;
        j0 = floor(ry_hy) - s_edge;
        k0 = floor(rz_hz) - s_edge;
        
        for (v = 0; v < p ; v++)
        {
            tx = rx_hx - i0 - v;
            ty = ry_hy - j0 - v;
            tz = rz_hz - k0 - v;
            phix[v] = base->region[v](tx);
            phiy[v] = base->region[v](ty);
            phiz[v] = base->region[v](tz);
        }
        
        for (i = 0; i < p ; i++)
        {
            for (j = 0 ; j < p ; j++)
            {
                for (k = 0; k < p ; k++)
                {
                    g = phix[i] * phiy[j] * phiz[k] * mass;
                    /* printf("%d %d %d --> %f\n",i+i0,j+j0,k+k0,g); */
                    gold = gridmass->getElement(gridmass,i+i0,j+j0,k+k0);
                    gridmass->setElement(gridmass,i+i0,j+j0,k+k0,g+gold);
                }
            }
        }
        curr = curr->next;
    }
}

void msm4g_grid_print(AbstractGrid *grid)
{
    int nx,ny,nz;
    int i,j,k;
    nx = grid->nx;
    ny = grid->ny;
    nz = grid->nz;
    for (k=0;k<nz;k++)
    {
        printf("k:%-2d\n",k);
        printf("%4s ","----");
        for (i=0; i<nx ; i++)
        {
            printf("i:%-7d ",i);
        }
        printf("\n");
        for (j=ny-1;j >= 0 ; j--)
        {
            printf("j:%-2d ",j);
            for (i=0; i<nx ; i++)
            {
                printf("%-9.5f ",grid->getElement(grid,i,j,k));
            }
            printf("\n");
        }
    }
}

void msm4g_grid_destroy(AbstractGrid **grid)
{
    (*grid)->destructor(grid);
    free(*grid);
    *grid = NULL;
}

AbstractGrid *msm4g_grid_dense_new(int nx, int ny, int nz,double h)
{
    AbstractGrid *grid;
    DenseGrid *densegrid;
    
    densegrid = malloc(sizeof(DenseGrid));
    grid = (AbstractGrid *)densegrid;
    grid->h  = h;
    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    grid->constructor = msm4g_grid_dense_new;
    grid->destructor  = msm4g_grid_dense_destroy;
    grid->setElement  = msm4g_grid_dense_setElement;
    grid->getElement  = msm4g_grid_dense_getElement;
    grid->reset       = msm4g_grid_dense_reset;

    densegrid->data = calloc(nx*ny*nz,sizeof(double));
    
    return grid;
}

void msm4g_grid_dense_setElement(AbstractGrid *grid,int i,int j,int k,double value)
{
    DenseGrid *densegrid = (DenseGrid *)grid;
    int position;
    
    position =  k * grid->nx * grid->ny + j * grid->nx  + i;
    densegrid->data[position] = value;
}

double msm4g_grid_dense_getElement(AbstractGrid *grid,int i,int j,int k)
{
    DenseGrid *densegrid = (DenseGrid *) grid;
    int position;
    
    position =  k * grid->nx * grid->ny + j * grid->nx  + i;
    return densegrid->data[position];
}

void msm4g_grid_dense_reset(AbstractGrid *grid,double value)
{
    int i,size;
    DenseGrid *self;
    self = (DenseGrid *)grid;
    size = grid->nx * grid->ny * grid->nz ;
    for (i=0; i<size; i++)
    {
        self->data[i] = value;
    }
}

void msm4g_grid_dense_destroy(AbstractGrid **grid)
{
    DenseGrid *densegrid;
    densegrid = (DenseGrid *)(*grid);
    free(densegrid->data);
    densegrid->data = NULL;
}

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

double msm4g_force_short(LinkedList *binlist,double threshold, msm4g_smoothing_handler smoothing_function)
{
    Bin *bin;
    Bin *neighbor;
    LinkedListElement *currBin;
    LinkedListElement *neighborBin;
    double potentialTotal = 0.0;
    double potential;
    
    currBin=binlist->head;
    while (currBin != NULL)
    {
        bin = (Bin *)currBin->data;
        potential = msm4g_force_short_withinBin(bin->particles,threshold,smoothing_function);
        potentialTotal += potential;
        
        neighborBin = bin->neighbors->head;
        while (neighborBin != NULL)
        {
            neighbor = (Bin *)neighborBin->data;
            if (neighbor->cantorindex > bin->cantorindex )
            {
                potential = msm4g_force_short_betweenBin(bin->particles,neighbor->particles,threshold,smoothing_function);
                potentialTotal += potential;
            }
            neighborBin=neighborBin->next;
        }
        currBin = currBin->next;
    }
    return potentialTotal;
}

double msm4g_force_short_withinBin(LinkedList *particles, double threshold, msm4g_smoothing_handler smoothing_function)
{
    Particle *particleI;
    Particle *particleJ;
    LinkedListElement *currI;
    LinkedListElement *currJ;
    double potentialTotal = 0.0;
    double potential;
    
    currI = particles->head;
    while (currI->next != NULL)
    {
        particleI = (Particle *)currI->data;
        currJ = currI->next;
        while (currJ != NULL)
        {
            particleJ = (Particle *)currJ->data;
            potential = msm4g_force_short_particlePair(particleI,particleJ,threshold,smoothing_function);
            potentialTotal += potential;
            currJ = currJ->next;
        }
        currI = currI->next;
    }
    return potentialTotal;
}

double msm4g_force_short_betweenBin(LinkedList *particlesI, LinkedList *particlesJ,double threshold, msm4g_smoothing_handler smoothing_function)
{
    LinkedListElement *currI, *currJ;
    Particle *particleI,*particleJ;
    double potential, potentialTotal = 0.0;
    currI = particlesI->head;
    while (currI != NULL)
    {
        particleI = (Particle *)currI->data;
        currJ = particlesJ->head;
        while (currJ != NULL)
        {
            particleJ = (Particle *)currJ->data;
            potential = msm4g_force_short_particlePair(particleI,particleJ,threshold,smoothing_function);
            potentialTotal += potential;
            currJ = currJ->next;
        }
        currI = currI->next;
    }
    return potentialTotal;
}

double msm4g_force_short_particlePair(Particle *particleI, Particle *particleJ, double a, msm4g_smoothing_handler gamma)
{
    D3Vector rij;
    double r,r2;
    double a2;
    double smoothedkernel,smoothedkernelderivative;
    double coeff;
    double potential = 0.0;
    
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
        coeff = particleI->m * particleJ->m * smoothedkernelderivative/r;
        /* f_j = f_j + coeff*rij */
        msm4g_d3vector_daxpy(&(particleJ->fshort) , &(particleJ->fshort), coeff, &rij);
        /* f_i = f_i - coeff*rij */
        msm4g_d3vector_daxpy(&(particleI->fshort) , &(particleI->fshort), coeff, &rij);
        potential = particleI->m * particleJ->m * smoothedkernel;
    }
    return potential;
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
        if (curr->data)
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
    msm4g_d3vector_set(&(box->width), 0.0,0.0,0.0);
    return box;
}

void msm4g_box_update(SimulationBox *box,LinkedList *list,double margin,double h,double p)
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
    
    /* enclosing the particles tightly */
    for (i=0;i<3;i++)
    {
        /* If the width in ith dimension is effectively zero,
         * we artificaially add some width (h) so that
         * a grid can be created.
         */
        if (fabs(max[i]-min[i]) < DBL_EPSILON )
        {
            box->location.value[i] = min[i]-h*0.50;
            box->width.value[i]    = h;
        } else
        {
            box->location.value[i] = min[i];
            box->width.value[i] = max[i]-min[i];
        }
    }
    
    /* enlarging the box <margin> percent. */
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
    
    /* make sure the widths are multiple of h */
    for (i=0;i<3;i++)
    {
        oldwidth = box->width.value[i];
        newwidth = fabs(floor(oldwidth/h)-oldwidth/h) < 0.0001 ? oldwidth : (floor(oldwidth/h)+1)*h;
        box->width.value[i] = newwidth;
        box->location.value[i] -= (newwidth-oldwidth)*0.5;
    }
    
    /* add p/2-1 boxes around the boundary */
    for (i=0;i<3;i++)
    {
        oldwidth = box->width.value[i];
        newwidth = oldwidth + (p-2)*h ;
        box->width.value[i] = newwidth;
        box->location.value[i] -= h*(p/2-1);
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
    msm4g_d3vector_set(&(particle->fshort), 0, 0, 0);
    msm4g_d3vector_set(&(particle->flong), 0, 0, 0);
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
            particle[j]->fshort.value[i] = (double)rand()/RAND_MAX;
            particle[j]->flong.value[i] = (double)rand()/RAND_MAX;
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
            particle = msm4g_particle_new(mass, r[0],r[1],r[2]);
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
        particle->fshort.value[i] = 0.0;
        particle->flong.value[i] = 0.0;
    }
    
    index++;
    return particle;
}

Particle *msm4g_particle_new(double mass,double x,double y,double z)
{
    Particle *particle;
    particle = msm4g_particle_empty();
    particle->m = mass;
    particle->r.value[0] = x;
    particle->r.value[1] = y;
    particle->r.value[2] = z;
    return particle;
}

void msm4g_particle_print(Particle *particle)
{
    printf("i:%d ",particle->index);
    printf("m:%8.2E ",particle->m);
    printf("r:%8.2E %8.2E %8.2E ",particle->r.value[0],particle->r.value[1],particle->r.value[2]);
    printf("v:%8.2E %8.2E %8.2E ",particle->v.value[0],particle->v.value[1],particle->v.value[2]);
    printf("fshort:%8.2E %8.2E %8.2E ",particle->fshort.value[0],particle->fshort.value[1],particle->fshort.value[2]);
    printf("flong:%8.2E %8.2E %8.2E ",particle->flong.value[0],particle->flong.value[1],particle->flong.value[2]);
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

