/** @file msm4g_lib.c
 * @brief The definitions of all core functions of MSM4G package.
 */
#include "msm4g_lib.h"
#include "msm4g_bases.h"

Simulation *msm4g_simulation_new(char *datafile,SimulationBox *box,Boolean periodic,int order,double abar,int mu) {
    Simulation *simulation;
    SimulationParameters *sp;
    SimulationOutput *output;
    Particle *particles;
    int N;
    
    simulation = (Simulation *)calloc(1,sizeof(Simulation));
    sp = (SimulationParameters *)calloc(1,sizeof(SimulationParameters));
    output = (SimulationOutput *)calloc(1,sizeof(SimulationOutput));
    particles = msm4g_particle_read(datafile,&N);
    sp->N = N;
    sp->abar = abar;
    sp->mu = mu;
    sp->periodic = periodic;
    sp->nu = order ;
    sp->L = floor(1 + 0.5 * log(sp->N)/log(8));
    double vol = box->wx * box->wy * box->wz ;
    sp->Mx = ceil(pow(2,1-sp->L)*box->wx * pow(sp->N/(vol),1/3.0)) ;
    sp->My = ceil(pow(2,1-sp->L)*box->wy * pow(sp->N/(vol),1/3.0)) ;
    sp->Mz = ceil(pow(2,1-sp->L)*box->wz * pow(sp->N/(vol),1/3.0)) ;

    sp->Mxmin = - (sp->Mx - 1) / 2;
    sp->Mxmax =    sp->Mx      / 2;
    sp->Mymin = - (sp->My - 1) / 2;
    sp->Mymax =    sp->My      / 2;
    sp->Mzmin = - (sp->Mz - 1) / 2;
    sp->Mzmax =    sp->Mz      / 2;

    sp->hx = box->wx / sp->Mx ;
    sp->hy = box->wy / sp->My ;
    sp->hz = box->wz / sp->Mz ;
    sp->h  = MSM4G_MAX3(sp->hx,sp->hy,sp->hz);
    sp->a = sp->abar * sp->h ;

    simulation->box = box;
    simulation->particles = particles ;
    simulation->parameters = sp;
    simulation->output = output;

    return simulation;
}

void msm4g_simulation_run(Simulation *simulation) {
    SimulationParameters *sp = simulation->parameters;
    SimulationBox *box = simulation->box;
    Particle *particles = simulation->particles;

    LinkedList *binlist = msm4g_bin_generate(box,particles,simulation->parameters->N,sp->a);
    
    msm4g_force_short(binlist, sp->a, simulation);
    msm4g_anterpolation(simulation);

    /* Calculate short-range potential energy */
    {
        double energy = 0.0;
        for (int i = 0 ; i < simulation->parameters->N ; i++) {
            energy += particles[i].potential_short_real * particles[i].m ;
        }
        simulation->output->potentialEnergyShortRange = energy;
    }


    msm4g_bin_destroy(binlist);
}

void msm4g_simulation_delete(Simulation *simulation) {
    if (simulation->grid != NULL)
        msm4g_grid_destroy(&(simulation->grid));
    msm4g_box_destroy(simulation->box);
    free(simulation->particles);
    free(simulation->parameters);
    free(simulation->output);
    free(simulation);
}

void msm4g_anterpolation(Simulation *simulation) //AbstractGrid *gridmass,SimulationBox *box,LinkedList *particles,const BaseFunction *base)
{
    int i0,j0,k0;
    int p = simulation->parameters->nu ;
    double x,y,z;
    double x0,y0,z0;
    double rx_hx,ry_hy,rz_hz;
    int s_edge;
    double h = simulation->parameters->h ;
    int mx = simulation->parameters->Mx ;
    int my = simulation->parameters->My ;
    int mz = simulation->parameters->Mz ;
    double tx,ty,tz;
    double phix[MAX_POLY_DEGREE], phiy[MAX_POLY_DEGREE], phiz[MAX_POLY_DEGREE];
    Particle *particles = simulation->particles;
    SimulationBox *box = simulation->box ;
    

    int nu = simulation->parameters->nu ;
    int N = simulation->parameters->N ;
    AbstractGrid *grid ;

    

    x0 = box->location.value[0];
    y0 = box->location.value[1];
    z0 = box->location.value[2];
    if (simulation->parameters->periodic == true) {
        grid = msm4g_grid_dense_new(mx,my,mz, h);
    } else {
        s_edge = p/2-1;
        grid = msm4g_grid_dense_new(mx+p-1,my+p-1,mz+p-1, h);
        for (int particleindex = 0 ; particleindex < N ; particleindex++) {
        
            Particle *particle = &(particles[particleindex]);
            x =  particle->r.value[0];
            y =  particle->r.value[1];
            z =  particle->r.value[2];
            double mass = particle->m;
            
            rx_hx = (x-x0) / h;
            ry_hy = (y-y0) / h;
            rz_hz = (z-z0) / h;
            
            i0 = floor(rx_hx) - s_edge;
            j0 = floor(ry_hy) - s_edge;
            k0 = floor(rz_hz) - s_edge;
            
            for (int v = 0; v < nu ; v++)
            {
                tx = rx_hx - i0 - v;
                ty = ry_hy - j0 - v;
                tz = rz_hz - k0 - v;
                phix[v] = msm4g_bases_bspline(nu,tx+nu/2) ;
                phiy[v] = msm4g_bases_bspline(nu,ty+nu/2) ;
                phiz[v] = msm4g_bases_bspline(nu,tz+nu/2) ;
            }
            
            for (int i = 0; i < nu ; i++)
            {
                for (int j = 0 ; j < nu ; j++)
                {
                    for (int k = 0; k < nu ; k++)
                    {
                        double g = phix[i] * phiy[j] * phiz[k] * mass;
                        double gold = grid->getElement(grid,i+i0,j+j0,k+k0);
                        grid->setElement(grid,i+i0,j+j0,k+k0,g+gold);
                    }
                }
            }
        }
    }
    simulation->grid = grid ;
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
    if (grid != NULL) {
        (*grid)->destructor(grid);
        free(*grid);
    }
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

double msm4g_smoothing_gama(double rho,int nu)
{
    double out;
    if (rho >= 1.0)
        out = 1.0 / rho;
    else {
        double rho2m1 = rho * rho - 1.0;
        out = 1.0;
        for (int k = nu - 1; k >= 1; k--) {
            out = 1.0 + ((0.5 - k) / k) * rho2m1 * out;
        }
    }
    return out;
}

double msm4g_smoothing_gamaprime(double rho,int nu)
{
    double outprime;
    double rho2 = rho * rho;
    if (rho >= 1.0) {
        outprime = -1.0 / rho2;
    } else {
        double rho2m1 = rho * rho - 1.0;
        outprime = nu-1;
        for (int k = nu - 1; k >= 2; k--) {
            outprime = k - 1 + (0.5-k)/k * rho2m1 * outprime;
        }
        outprime = - rho * outprime ;
    }
    return outprime ;
}

double msm4g_force_short(LinkedList *binlist,double threshold, Simulation *simulation)
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
        potential = msm4g_force_short_withinBin(bin->particles,threshold,simulation);
        potentialTotal += potential;
        
        neighborBin = bin->neighbors->head;
        while (neighborBin != NULL)
        {
            neighbor = (Bin *)neighborBin->data;
            if (neighbor->cantorindex > bin->cantorindex )
            {
                potential = msm4g_force_short_betweenBin(bin->particles,neighbor->particles,threshold,simulation);
                potentialTotal += potential;
            }
            neighborBin=neighborBin->next;
        }
        currBin = currBin->next;
    }
    return potentialTotal;
}

double msm4g_force_short_withinBin(LinkedList *particles, double threshold, Simulation *simulation)
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
            potential = msm4g_force_short_particlePair(particleI,particleJ,threshold,simulation);
            potentialTotal += potential;
            currJ = currJ->next;
        }
        currI = currI->next;
    }
    return potentialTotal;
}

double msm4g_force_short_betweenBin(LinkedList *particlesI, LinkedList *particlesJ,double threshold, Simulation *simulation)
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
            potential = msm4g_force_short_particlePair(particleI,particleJ,threshold,simulation);
            potentialTotal += potential;
            currJ = currJ->next;
        }
        currI = currI->next;
    }
    return potentialTotal;
}

double msm4g_force_short_particlePair(Particle *particleI, Particle *particleJ, double a, Simulation *simulation)
{
    D3Vector rij;
    double r,r2;
    double a2 = a * a;
    double smoothedkernel,smoothedkernelderivative;
    double coeff;
    double potential = 0.0;
    int order = simulation->parameters->nu;
    Boolean periodic = simulation->parameters->periodic;
    if (periodic == true) {
        double a = simulation->parameters->a ;
        double Ax = simulation->box->wx;
        double Ay = simulation->box->wy;
        double Az = simulation->box->wz;
        double rix = particleI->r.value[0] ;
        double riy = particleI->r.value[1] ;
        double riz = particleI->r.value[2] ;
        double rjx = particleJ->r.value[0] ;
        double rjy = particleJ->r.value[1] ;
        double rjz = particleJ->r.value[2] ;
        int pxmin = (rix - rjx - a)/Ax ;
        int pxmax = (rix - rjx + a)/Ax ;
        int pymin = (riy - rjy - a)/Ay ;
        int pymax = (riy - rjy + a)/Ay ;
        int pzmin = (riz - rjz - a)/Az ;
        int pzmax = (riz - rjz + a)/Az ;
        for (int px = pxmin ; px <= pxmax ; px++) {
            double r2i =  (rix-rjx-Ax*px)*(rix-rjx-Ax*px) ;
            for (int py = pymin ; py <= pymax ; py++) {
                double r2j = r2i + (riy-rjy-Ay*py)*(riy-rjy-Ay*py);
                for (int pz = pzmin ; pz <= pzmax ; pz++) {
                    double r2 = r2j + (riz-rjz-Az*pz)*(riz-rjz-Az*pz);
                    double r = sqrt(r2);
                    double g0 =  1.0/r -  msm4g_smoothing_gama(r/a,order)/a;
                    particleI->potential_short_real += 0.5 * particleJ->m  * g0;
                    particleJ->potential_short_real += 0.5 * particleI->m  * g0;
                    double rx = rix - rjx - Ax * px ;
                    double ry = riy - rjy - Ay * py ;
                    double rz = riz - rjz - Az * pz ;
                    double gamaprime = (-1./r2 - (1./a2) * msm4g_smoothing_gamaprime(r/a, order))/r ;
                    particleI->acc_short[0] +=  particleJ->m * gamaprime * rx ;
                    particleI->acc_short[1] +=  particleJ->m * gamaprime * ry ;
                    particleI->acc_short[2] +=  particleJ->m * gamaprime * rz ;
                    particleJ->acc_short[0] -=  particleI->m * gamaprime * rx ;
                    particleJ->acc_short[1] -=  particleI->m * gamaprime * ry ;
                    particleJ->acc_short[2] -=  particleI->m * gamaprime * rz ;
                }
            }
        }

    } else {
        /* rij = rj - ri */
        msm4g_d3vector_daxpy(&rij, &(particleJ->r), -1.0, &(particleI->r));
        /* r^2 = |rij|^2 */
        r2=msm4g_d3vector_normsquare(&rij);
        r=sqrt(r2);
        /* If only the particles are closer than cut-off distance */
        if (r2 < a2)
        {
            smoothedkernel = 1.0/r -  msm4g_smoothing_gama(r/a,order)/a;
            smoothedkernelderivative = 1.0/r2 - msm4g_smoothing_gamaprime(r/a,order)/a2;
            coeff = particleI->m * particleJ->m * smoothedkernelderivative/r;
            /* f_j = f_j + coeff*rij */
            msm4g_d3vector_daxpy(&(particleJ->fshort) , &(particleJ->fshort), coeff, &rij);
            /* f_i = f_i - coeff*rij */
            msm4g_d3vector_daxpy(&(particleI->fshort) , &(particleI->fshort), coeff, &rij);
            potential = particleI->m * particleJ->m * smoothedkernel;
            particleI->potential_short_real += 0.5 * particleJ->m  * smoothedkernel;
            particleJ->potential_short_real += 0.5 * particleI->m  * smoothedkernel;
        }
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

SimulationBox *msm4g_box_newCube(double location, double width)
{
    SimulationBox *box = msm4g_box_new();
    box->location.value[0] = location;
    box->location.value[1] = location;
    box->location.value[2] = location;
    box->width.value[0] = width;
    box->width.value[1] = width;
    box->width.value[2] = width;
    box->x  = box->y  = box->z  = location;
    box->wx = box->wy = box->wz = width;
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

LinkedList *msm4g_bin_generate(SimulationBox *box,Particle *particles,int n,double binwidth)
{
    LinkedList *binlist;
    Bin *bin;
    I3Vector binindex;

    binlist = msm4g_linkedlist_new();
    
    for (int particleindex = 0 ; particleindex < n ; particleindex++ )
    {
        Particle *particle = &(particles[particleindex]);
        msm4g_i3vector_set(&binindex, -1, -1, -1);
        for (int i=0;i<3;i++)
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

Particle *msm4g_particle_read(const char *filename,int *numberofparticles)
{
    Particle *particles;
    FILE *fp;
    int n;
    double mass;
    double r[3];
    double v[3];
    
    fp = fopen(filename,"r");
    if (fp == NULL)
    {
        return NULL;
    }
    fscanf(fp,"%d",&n);
    *numberofparticles = n;
    particles = (Particle *)calloc(n,sizeof(Particle));
    for (int iparticle = 0 ; iparticle < n ; iparticle++)
    {
        int ret = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&mass,&r[0],&r[1],&r[2],&v[0],&v[1],&v[2]);
        if (ret == 7)
        {
            particles[iparticle].m = mass;
            particles[iparticle].r.value[0] = r[0];
            particles[iparticle].r.value[1] = r[1];
            particles[iparticle].r.value[2] = r[2];
            particles[iparticle].v.value[0] = v[0];
            particles[iparticle].v.value[1] = v[1];
            particles[iparticle].v.value[2] = v[2];

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
        particle->acc_short[i] = 0.0;
    }

    particle->potential_short_real = 0.0;
    
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
    /* printf("i:%d ",particle->index);
    printf("m:%8.2E ",particle->m);
    printf("r:%8.2E %8.2E %8.2E ",particle->r.value[0],particle->r.value[1],particle->r.value[2]);
    printf("v:%8.2E %8.2E %8.2E ",particle->v.value[0],particle->v.value[1],particle->v.value[2]); */
    printf("%25.16E %25.16E %25.16E\n",particle->acc_short[0],particle->acc_short[1],particle->acc_short[2]);
    /* printf("flong:%8.2E %8.2E %8.2E ",particle->flong.value[0],particle->flong.value[1],particle->flong.value[2]); */
    /* printf("\n");*/
    
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

double *msm4g_util_omegaprime(int mu,int p) {
    int iteration;
    if (p != 4 && p != 6) {
        fprintf(stderr, "p can be 4 or 6 only\n");
    }
    double *wprime = (double *)calloc(mu+p/2+1, sizeof(double));
    
    double P[100], gold[100], g[100], gprime[100], ksiprime[100], psiprime[100];
    double Aprime[100];
    // AppendixC --> Cpoly outputs only for p=4 and p=6
    double psiprime4[2] = {  6.4814814814814811e-02,  9.2592592592592587e-03 };
    double psiprime6[4] = { -3.0804398148148150e-02, -8.9800347222222226e-03,
        -6.0590277777777780e-04, -1.1863425925925925e-05 };
    double firstpart4[3] = {-3.3333333333333331e-01, 1.6666666666666665e+00,
        -3.3333333333333331e-01 };
    double firstpart6[5] = { 1.7083333333333334e-01, -1.1833333333333333e+00,
        3.0249999999999999e+00, -1.1833333333333333e+00,
        1.7083333333333334e-01 };
    double deltastencil[100], stencil[100];
    //printf("mu = %d\n",mu);
    for (int i = 0; i < 100; i++) {
        P[i] = 0.0;
        gold[i] = 0.0;
        g[i] = 0.0;
        gprime[i] = 0.0;
        ksiprime[i] = 0.0;
        psiprime[i] = 0.0;
        deltastencil[i] = 0.0;
        stencil[i] = 0.0;
    }
    
    for (int i = 0; i < p / 2; i++) {
        P[i] = msm4g_bases_bspline(p, i + p / 2);
        gold[i] = 0.0;
        g[i] = 0.0;
        //printf("P[%d] %f\n",i,P[i]);
    }
    gold[0] = 1;
    iteration = 0;
    while (1) {
        double value = 0.0;
        for (int k = 1; k < p / 2; k++)
            value += gold[k] * gold[k];
        value = P[0] - value;
        if (value < 0) {
            fprintf(stderr, "something wrong %d\n", iteration);
            return wprime;
        }
        g[0] = sqrt(value);
        for (int j = 1; j <= p / 2 - 1; j++) {
            double sumtemp = 0.0;
            for (int i = j + 1; i <= p / 2 - 1; i++) {
                sumtemp += gold[i - j] * gold[i];
            }
            g[j] = (P[j] - sumtemp) / g[0];
        }
        double relerr = msm4g_util_diffnorm(g, gold, p / 2) / msm4g_util_norm(g, p / 2);
        if (relerr < 1e-16 || iteration > 1000)
            break;
        for (int k = 0; k < p / 2; k++)
            gold[k] = g[k];
        iteration++;
    }
    //for (int k=0;k<p/2;k++)
    //printf("g[%d] %f\n",k,g[k]);
    
    // Calculating gprime G2=GG
    for (int k = 0; k < p - 1; k++)
        gprime[k] = 0.0;
    for (int i = 0; i <= p - 2; i++) {
        int start = 0;
        if (i > p / 2 - 1)
            start = i - (p / 2 - 1);
        int stop = fmin(i, p / 2 - 1);
        for (int j = start; j <= stop; j++)
            gprime[i] += g[j] * g[i - j];
    }
    //for (int k=0;k <p-1;k++)
    //  printf("gprime[%d] %25.16e\n",k,gprime[k]);
    for (int i = 0; i < p - 2; i++) {
        if (p == 4)
            psiprime[i] = psiprime4[i];
        else if (p == 6)
            psiprime[i] = psiprime6[i];
    }
    
    for (int i = 0; i < p - 1; i++)
        ksiprime[i] = 0.0;
    for (int i = p - 3; i >= 0; i--) {
        double sumprime = 0.0;
        for (int j = i + 1; j <= p - 3; j++) {
            sumprime += ksiprime[j] * gprime[j - i];
        }
        ksiprime[i] = (psiprime[i] - sumprime) / gprime[0];
    }
    //for (int i = 0 ; i < p-1 ; i++)
    //  printf("ksiprime[%d] = %25.16e\n",i,ksiprime[i]);
    for (int i = 0; i < p - 1; i++)
        for (int j = 0; j < p - 1; j++)
            Aprime[i * (p - 1) + j] = 0.0;
    for (int i = 0; i <= p - 2; i++) {
        for (int j = 0; j <= p - 2; j++) {
            if (i >= j)
                Aprime[i * (p - 1) + j] = gprime[i - j];
        }
    }
    for (int j = 1; j <= p - 2; j++) {
        for (int i = 0; i <= p - 2 - j; i++) {
            Aprime[i * (p - 1) + j] += gprime[i + j];
        }
    }
    /*
     for (int i=0;i<p-1;i++) {
     for (int j=0;j<p-1;j++)
     printf("%25.16f ",Aprime[i*(p-1)+j]);
     printf(" | %8.5f\n",ksiprime[i]);
     }*/
    
    msm4g_util_gausssolver(p - 1, Aprime, ksiprime);
    /*
     printf("after\n");
     for (int i=0;i<p-1;i++) {
     for (int j=0;j<p-1;j++)
     printf("%25.16f ",Aprime[i*(p-1)+j]);
     printf(" | %8.5f\n",ksiprime[i]);
     } */
    
    double c[100];
    for (int i = 0; i < 100; i++)
        c[i] = 0;
    for (int i = 0; i < p - 1; i++)
        c[i] = ksiprime[i];
    for (int n = p - 1; n <= 99; n++) {
        double tmpsum = 0.0;
        for (int k = 1; k <= p - 2; k++) {
            tmpsum += c[n - k] * gprime[k];
        }
        c[n] = -tmpsum / gprime[0];
    }
    /* for (int i = 0 ; i < 100 ; i++)
     printf("c[3d] %18.15f\n",c[i]); */
    
    for (int i = 0; i <= p - 1; i++) {
        if (p == 4)
            deltastencil[i] = firstpart4[i];
        else if (p == 6)
            deltastencil[i] = firstpart6[i];
    }
    //for (int i = mu + p/2  ; i <= 2*mu+p ;i++ )
    //    printf("stencil[%3d] : %18.15f\n",i-mu-p/2,stencil[i]);
    
    for (int i = 1; i <= p - 1; i++) {
        double value = deltastencil[i - 1];
        int location = i - p / 2 + mu + p / 2;
        stencil[location] += value;
    }
    /* for (int i = 0 ; i < 2*mu+p+1 ; i++)
     printf("stencil[%3d] : %18.16f\n",i,stencil[i]);*/
    
    double delta2stencil4[5] = { 1, -4, 6, -4, 1 };
    double delta2stencil6[7] = { 1, -6, 15, -20, 15, -6, 1 };
    double delta2stencil[50];
    for (int i = 0; i < p + 1; i++) {
        if (p == 4)
            delta2stencil[i] = delta2stencil4[i];
        else if (p == 6)
            delta2stencil[i] = delta2stencil6[i];
    }
    
    for (int m = -mu; m <= mu; m++) {
        for (int i = 1; i <= p + 1; i++) {
            double value = delta2stencil[i - 1] * c[abs(m)];
            int location = i - p / 2 - 1 + m + mu + p / 2;
            stencil[location] += value;
        }
    }
    for (int i = mu + p / 2; i <= 2 * mu + p; i++) {
        //printf("stencil[%3d] : %18.15f\n",i-mu-p/2,stencil[i]);
        wprime[i - mu - p / 2] = stencil[i];
    }
    //for (int i=0;i<20;i++)
    // printf("wprime[%3d] : %18.15f\n",i,wprime[i]);
    return wprime;
}

double msm4g_util_norm(double x[], int n) {
    double out = 0.0;
    for (int i = 0; i < n; i++)
        out += x[i] * x[i];
    out = sqrt(out);
    return out;
}

double msm4g_util_diffnorm(double x[], double y[], int n) {
    double out = 0.0;
    for (int i = 0; i < n; i++)
        out += (x[i] - y[i]) * (x[i] - y[i]);
    return sqrt(out);
}

void msm4g_util_gausssolver(int n, double *A, double *y) {
    /* Forward substition */
    for (int j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            double multiplier = -A[i * n + j] / A[j * n + j];
            for (int k = 0; k < n; k++) {
                A[i * n + k] += multiplier * A[j * n + k];
            }
            y[i] += multiplier * y[j];
        }
    }
    /* Backward substition */
    for (int j = n - 1; j >= 0; j--) {
        y[j] = y[j] / A[j * n + j];
        A[j * n + j] = 1.0;
        for (int i = j - 1; i >= 0; i--) {
            y[i] = y[i] - A[i * n + j] * y[j];
            A[i * n + j] = 0.0;
        }
    }
}
