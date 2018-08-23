/** @file msm4g_tests.c
 * @brief The test functions to check the functionality of the MSM4G package.
 *
 * This is probably the best place one can learn how to use the package.
 * These test functions are used to check the outcomes of the functions
 * by comparing with expected results.
 */
#include "msm4g_lib.h"
#include "msm4g_tests.h"

void msm4g_unit_test_all()
{
    Boolean status;
    int index, numberoftests,failed;
    typedef Boolean (testFunctionType)();
    testFunctionType *testFunction;
    LinkedList *list = msm4g_linkedlist_new();
    
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_1);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_2);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_3);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_4);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_5);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_6);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_7);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_8);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_9);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_10);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_11);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_12);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_13);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_14);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_15);

    numberoftests = msm4g_linkedlist_size(list);

    failed=0;
    for (index = 0; index < numberoftests ; index++)
    {
        testFunction = (testFunctionType *)msm4g_linkedlist_get(list, index);
        status = testFunction();
        if (status == false)
        {
            failed++;
            fprintf(stderr,"Unit test %3d: %s\n",index+1,"failed");
        }
    }
    printf("%d of %d passed, %d failed\n",numberoftests-failed,numberoftests,failed);
    
    msm4g_linkedlist_destroy(list);

}

Boolean msm4g_unit_test_1()
{
    LinkedList *list;
    Particle x[10];
    int i;
    
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        msm4g_particle_reset(x+i);
        msm4g_linkedlist_add(list,&x[i]);
    }
    if (msm4g_linkedlist_size(list) != 10)
    {
        msm4g_linkedlist_destroy(list);
        return false;
    }
    msm4g_linkedlist_destroy(list);
    return true;
}

Boolean msm4g_unit_test_2()
{
    LinkedList *list;
    LinkedListElement *curr;
    Particle x[10];
    Particle *y;
    Boolean status ;
    int i;
    
    status = true;
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        x[i].m = i;
        msm4g_linkedlist_add(list,&x[i]);
    }
    
    /**< Iterating from tail to head */
    i=9;
    curr=list->tail;
    while (curr != NULL)
    {
        y =  (Particle *)curr->data;
        if (fabs(i - (y->m)) > DBL_EPSILON)
        {
            status = false;
            break;
        }
        curr = curr->prev;
        i--;
    }
    msm4g_linkedlist_destroy(list);
    
    return status;
}

Boolean msm4g_unit_test_3()
{
    LinkedList *list;
    LinkedListElement *curr;
    Particle x[10];
    Particle *y;
    Boolean status ;
    int i;
    
    status = true;
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        x[i].m = i;
        msm4g_linkedlist_add(list,&x[i]);
    }
    
    i=0;
    curr=list->head;
    while (curr != NULL)
    {
        y =  (Particle *)curr->data;
        if (fabs(i - (y->m)) > DBL_EPSILON)
        {
            status = false;
            break;
        }
        curr = curr->next;
        i++;
    }
    msm4g_linkedlist_destroy(list);
    
    return status;
}


Boolean msm4g_unit_test_4()
{
    const int DIM = 3;
    const int N = 3; /* Number of particles in eight.ini */
    int dummy;
    Particle particles[3];
    Particle *particlesInFile;
    int i,j;
    
    for (i=0;i<N;i++)
    {
        msm4g_particle_reset(&particles[i]);
        particles[i].m = 1.0;
    }
    particles[0].r.value[0] =   0.97000436 ;
    particles[0].r.value[1] =  -0.24308753 ;
    particles[0].v.value[0] =   0.466203685 ;
    particles[0].v.value[1] =   0.43236573  ;
    particles[1].r.value[0] =  -0.97000436;
    particles[1].r.value[1] =   0.24308753;
    particles[1].v.value[0] =   0.466203685;
    particles[1].v.value[1] =   0.43236573;
    particles[2].v.value[0] =  -0.93240737;
    particles[2].v.value[1] =  -0.86473146;
    
    particlesInFile=msm4g_particle_read("data/eight.ini",&dummy);
   
    
    for (i=0;i<N;i++) /* for each particle */
    {
        if (fabs(particles[i].m - particlesInFile[i].m)>DBL_EPSILON ) {
            return false;
        }
        for (j=0;j<DIM;j++) /* for each dimension */
        {
            if (fabs(particles[i].r.value[j]-particlesInFile[i].r.value[j]) > DBL_EPSILON) {
                return false;
            }
            if (fabs(particles[i].v.value[j]-particlesInFile[i].v.value[j]) > DBL_EPSILON) {
                return false;
            }
        }
    }
    
    /* Deallocate the particles */
    free(particlesInFile);
    return true;
}

Boolean msm4g_unit_test_5()
{
    int i,n,p;
    Boolean status ;
    LinkedList *list ;
    Particle **particles;
    SimulationBox *box;
    D3Vector locationError;
    D3Vector locationExpected ;
    D3Vector widthExpected;
    D3Vector widthError;
    double h;
    
    
    msm4g_d3vector_set(&locationExpected, -1.5, -3.0, -4.5);
    msm4g_d3vector_set(&widthExpected,     3.0,  6.0,  9.0);
    
    n=10;
    status = true;
    particles = msm4g_particle_rand(n);
    
    list = msm4g_linkedlist_new();
    for (i=0;i<n;i++)
    {
        msm4g_linkedlist_add(list, particles[i]);
    }
    
    particles[5]->r.value[0] = -1.0;
    particles[5]->r.value[1] = -2.0;
    particles[5]->r.value[2] = -3.0;
    
    particles[7]->r.value[0] =  1.0;
    particles[7]->r.value[1] =  2.0;
    particles[7]->r.value[2] =  3.0;
    
    box = msm4g_box_new();
    h=1.0; /* Lattice spacing in the finest level */
    p=2; /* degree of base polynomials. if p=1 there is no need for padding around the boundary */
    msm4g_box_update(box, list, 0.5,h,p);
    
    msm4g_d3vector_daxpy(&locationError, &locationExpected,-1.0,&(box->location));
    msm4g_d3vector_daxpy(&widthError,    &widthExpected,   -1.0,&(box->width));
    
    if (msm4g_d3vector_norm(&locationError) > DBL_EPSILON ||
        msm4g_d3vector_norm(&widthError)    > DBL_EPSILON )
    {
        status = false;
    }
    
    msm4g_linkedlist_destroy(list);
    msm4g_particle_destroyarray(particles, n);
    msm4g_box_destroy(box);
    
    return status;
}

Boolean msm4g_unit_test_6()
{
    Boolean status = true;
    Particle *particlelist;
    LinkedList *binlist;
    Bin *bin;
    SimulationBox box;
    double binwidth = 10;
    int i,n;
    
    particlelist = msm4g_particle_read("data/bintest.ini",&n);
    for (i=0;i<3;i++)
    {
        box.location.value[i] = 0.0;
        box.width.value[i]    = 30.0;
    }
    binlist=msm4g_bin_generate(&box,particlelist,6,binwidth);
    
    bin = (Bin *)msm4g_linkedlist_get(binlist, 0);
    if (bin->cantorindex != 8) return false;
    /* It should has only one neighbor */
    if (msm4g_linkedlist_size(bin->neighbors) != 1) return false;
    /* Its neighbor's cantor index should be 17 */
    if (((Bin *)msm4g_linkedlist_get(bin->neighbors,0))->cantorindex != 17) return false;
    /* It should contain two particles */
    if (msm4g_linkedlist_size(bin->particles) != 2) return false;
    
    bin = (Bin *)msm4g_linkedlist_get(binlist, 1);
    /* This bin should has a neighbor */
    if (msm4g_linkedlist_size(bin->neighbors) != 0) return false;
    /* Its cantor index should be 25 */
    if (bin->cantorindex != 25) return false;
    /* It should contain only one particles */
    if (msm4g_linkedlist_size(bin->particles) != 1) return false;
    
    bin = (Bin *)msm4g_linkedlist_get(binlist, 2);
    /* It should has only one neighbor */
    if (msm4g_linkedlist_size(bin->neighbors) != 1) return false;
    /* Its neighbor's cantor index should be 8 */
    if (((Bin *)msm4g_linkedlist_get(bin->neighbors,0))->cantorindex != 8) return false;
    /* It should contain three particles */
    if (msm4g_linkedlist_size(bin->particles) != 3) return false;
    
    msm4g_bin_destroy(binlist);
    free(particlelist);
    return status;
}

Boolean msm4g_unit_test_7()
{
    Boolean status = true;
    int i,j;
    int inp[6] = {0, 1, 2, 3,  4, 5  };
    int out[6] = {1, 1, 2, 6, 24, 120};
    int x[4][4];
    int xexpected[4][4] = {{0,1,3,6},{2,4,7,11},{5,8,12,17},{9,13,18,24}};
    int vector2d[2];
    
    
    
    /* Checking factorial function */
    for (i=0; i<6; i++)
        if (msm4g_math_factorial(inp[i]) != out[i])
            return false;
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            vector2d[0] = i; vector2d[1]=j;
            x[i][j] = msm4g_math_cantor(vector2d,2);
            if (x[i][j] != xexpected[i][j]) return false;
        }
    }
    return status;
}

Boolean msm4g_unit_test_8()
{
    Boolean status = true;
    int i;
    int n = 3;
    
    for (i=0; i<n; i++)
    {
        /* \gamma(1.0)  = 1.0 for all gamma function. */
        if (fabs(msm4g_smoothing_gama(1,i+1)-1.0) > DBL_EPSILON) return false;
        /* \gamma'(1.0) = -1.0 for all gamma function. */
        if (fabs(msm4g_smoothing_gamaprime(1,i+1)+1.0) > DBL_EPSILON) return false;
        /* \gamma'(0.0) = 0.0 for all gamma function. */
        if (fabs(msm4g_smoothing_gamaprime(0,i+1)-0.0) > DBL_EPSILON) return false;

    }
    /* \gamma(0.0) = 3/2 for C1 gamma */
    if (fabs(msm4g_smoothing_gama(0,2)-3.0/2.0)   > DBL_EPSILON) return false;
    /* \gamma(0.0) = 15/8 for C2 gamma */
    if (fabs(msm4g_smoothing_gama(0,3)-15.0/8.0)  > DBL_EPSILON) return false;
    /* \gamma(0.0) = 35/16 for C3 gamma */
    if (fabs(msm4g_smoothing_gama(0,4)-35.0/16.0) > DBL_EPSILON) return false;

    return status;
}

Boolean msm4g_unit_test_9()
{
    Boolean teststatus = true;
    Boolean periodic = false;
    int order = 2 ;
    double abar = 1;
    int mu = 2;
    char *datafile = "data/bintest.ini";
    double boxlocation = 0.0;
    double boxwidth = 30.0;
    SimulationBox *box = msm4g_box_newCube(boxlocation,boxwidth);

    Simulation *simulation = msm4g_simulation_new(datafile,box,periodic,order,abar,mu);

    simulation->parameters->a = 10;

    msm4g_simulation_run(simulation);

    double potential = simulation->output->potentialEnergyShortRange ;
    double potentialExpected =  1127.0/6000 + (5000*sqrt(2.0)-1776.0)/4000.0 + (2000*sqrt(82.0)-17876.0)/164000.0;

    if (fabs(potentialExpected-potential)/potential  > 1E-15) {
        teststatus = false;
    }

    msm4g_simulation_delete(simulation);
    return teststatus;
}

Boolean msm4g_unit_test_10()
{
    Boolean status = true;
    AbstractGrid *grid;
    
    double h = 1.0;
    int nx = 3;
    int ny = 3;
    int nz = 3;
        
    grid = msm4g_grid_dense_new(nx,ny,nz,h);
    /* Check if it could allocated the object */
    if (grid == NULL) return false;
    
    grid->reset(grid,0.0);
    
    grid->setElement(grid,0,1,2,10.0); 
    
    msm4g_grid_destroy(&grid);
    /* The grid should point to NULL after destruction */
    if (grid != NULL) return false;
    
    return status;
}

Boolean msm4g_unit_test_11()
{
    Boolean teststatus = true;
    AbstractGrid *grid ;
    SimulationBox *box;
    double h=1.5;
    char *datafile = "data/singleton.ini";
    Boolean periodic = false;
    double phix[4] = {1, 23, 23, 1};
    double abar = 4;
    int mu = 4;

    box = msm4g_box_newCube(-2.25,4.5);
    Simulation *simulation = msm4g_simulation_new(datafile,box,periodic,4,abar,mu);
    simulation->parameters->Mx = 1;
    simulation->parameters->My = 1;
    simulation->parameters->Mz = 1;
    simulation->parameters->h  = h ;

    msm4g_simulation_run(simulation);
    
    grid = simulation->grid;
    for (int i = 0 ; i < grid->nx ; i++)
        for (int j = 0 ; j < grid->ny ; j++)
            for (int k = 0 ; k < grid->nz ; k++)
                if (fabs(grid->getElement(grid,i,j,k)-phix[i]*phix[j]*phix[k]) > FLT_EPSILON) {
                    teststatus = false;
                    break;
                }
    

    msm4g_simulation_delete(simulation);
    return teststatus;
}

Boolean msm4g_unit_test_12()
{
    /* Expected values are from Table I of doi:10.1063/1.4943868 */
    char *omegaprimeCubicExpected[13] = {"3.464","-1.732", "0.679","-0.240","0.080","-0.026","0.008","-0.002","0.001","-0.000","0.000","-0.000","0.000"};
    char *omegaprimeQuinticExpected[13] = {"12.379","-9.377", "5.809","-3.266","1.735","-0.889","0.444","-0.217","0.105","-0.050","0.024","-0.011","0.005"};
    double *omegaprimeCubicCalculated = msm4g_util_omegaprime(20, 4);
    double *omegaprimeQuinticCalculated = msm4g_util_omegaprime(20, 6);
    for (int i = 0 ; i < 13 ; i++) {
        char *expectedCubic   = omegaprimeCubicExpected[i];
        char *expectedQuintic = omegaprimeQuinticExpected[i];
        char calculated[7]  ;
        sprintf(&(calculated[0]),"%5.3f",omegaprimeCubicCalculated[i]);
        if (strcmp(expectedCubic, calculated) != 0) return false ;
        sprintf(&(calculated[0]),"%5.3f",omegaprimeQuinticCalculated[i]);
        if (strcmp(expectedQuintic, calculated) != 0) return false ;
    }
    free(omegaprimeCubicCalculated);
    free(omegaprimeQuinticCalculated);
    return true;
}

Boolean msm4g_unit_test_13()
{
    const int ntest = 9;
    const int cubic = 4;
    const int quintic = 6;
    const double mathematicaCubic[9]        = {0 ,1,   8,   23,   32,   23,   8,  1,0}; /* divided by    48 */
    const double mathematicaCubicPrime[9]   = {0 ,1,   4,    5,    0,   -5,  -4, -1,0}; /* divided by     8 */
    const double mathematicaQuintic[9]      = {0,81,2528,13438,22528,13438,2528, 81,0}; /* divided by 40960 */
    const double mathematicaQuinticPrime[9] = {0,27, 400,  942,    0, -942,-400,-27,0}; /* divided by  2048 */

    for (int i = 0 ; i < ntest ; i++) {
      double expectedCubic        = mathematicaCubic[i]       /   48.0 ;
      double expectedCubicPrime   = mathematicaCubicPrime[i]  /    8.0 ;
      double expectedQuintic      = mathematicaQuintic[i]     /40960.0 ;
      double expectedQuinticPrime = mathematicaQuinticPrime[i]/ 2048.0 ;
      double calculatedCubic        = msm4g_bases_bspline     (cubic  ,cubic  *i/(double)(ntest-1));
      double calculatedCubicPrime   = msm4g_bases_bsplineprime(cubic  ,cubic  *i/(double)(ntest-1));
      double calculatedQuintic      = msm4g_bases_bspline     (quintic,quintic*i/(double)(ntest-1));
      double calculatedQuinticPrime = msm4g_bases_bsplineprime(quintic,quintic*i/(double)(ntest-1));
      if (fabs(expectedCubic       -calculatedCubic)        > DBL_EPSILON) return false;
      if (fabs(expectedCubicPrime  -calculatedCubicPrime)   > DBL_EPSILON) return false;
      if (fabs(expectedQuintic     -calculatedQuintic)      > DBL_EPSILON) return false;
      if (fabs(expectedQuinticPrime-calculatedQuinticPrime) > DBL_EPSILON) return false;
    }
    return true;
}

Boolean msm4g_unit_test_14()
{
    double expectedgama[5][4] =  {
            {  3/  2.0,     11/     8.0, 1, 2/3.0}, /* nu = 2 */
            { 15/  8.0,    203/   128.0, 1, 2/3.0}, /* nu = 3 */
            { 35/ 16.0,   1759/  1024.0, 1, 2/3.0}, /* nu = 4 */
            {315/128.0,  59123/ 32768.0, 1, 2/3.0}, /* nu = 5 */
            {693/256.0, 488293/262144.0, 1, 2/3.0}};/* nu = 6 */
    double expectedgamaprime[5][4] = {
            {0,     -1/    2.0, -1, -4/9.0},
            {0,    -17/   16.0, -1, -4/9.0},
            {0,   -407/  256.0, -1, -4/9.0},
            {0,  -4201/ 2048.0, -1, -4/9.0},
            {0,-159947/65536.0, -1, -4/9.0}};
    for (int nu = 2; nu <= 6 ; nu++) {
        for (int i = 0 ; i < 4 ; i++) {
            double rho = i / 2.0 ;
            double gama = msm4g_smoothing_gama(rho,nu);
            double gamaprime = msm4g_smoothing_gamaprime(rho,nu);
            if (fabs(gama     -expectedgama[nu-2][i])      > DBL_EPSILON) return false;
            if (fabs(gamaprime-expectedgamaprime[nu-2][i]) > DBL_EPSILON) return false;
        }
    }
    return true;
}

Boolean msm4g_unit_test_15()
{
    Boolean teststatus = true;
    Boolean periodic = true;
    int order = 4;
    int mu = 2;
    double abar = 4.0;
    double expected = 6.9843177228211383e-09;
    /* double expectedacc[8][3] = {
        {   1.7127138343151520e-05,   -7.9922393447316331e-06,    7.0828748251281019e-06},
        {  -3.1966280586136999e-05,   -1.0376534820290217e-04,   -6.0893167284709774e-05},
        {  -4.3875970537235183e-05,   -2.7824718658855566e-05,    3.0876280003449423e-05},
        {   7.3727156274532452e-05,    1.1860770584138265e-04,   -1.1113166484930614e-05},
        {   7.8059074802170574e-06,   -1.4019342701374309e-07,   -5.7782346576733895e-05},
        {   5.9620212096497241e-05,   -2.3746627592581602e-05,    3.1717687081769986e-05},
        {   1.4998005198020302e-05,    2.9271945320370181e-05,    6.7963357688293238e-06},
        {  -9.7436168269046449e-05,    1.5589476064332005e-05,    5.3315502667197417e-05}}; */
    SimulationBox *unitCube = msm4g_box_newCube(0,1);

    Simulation *simulation = msm4g_simulation_new("data/changaN8.ini",unitCube,periodic,order,abar,mu);

    msm4g_simulation_run(simulation);
    
    double relativeError = fabs(simulation->output->potentialEnergyShortRange-expected)/fabs(expected);
    if (relativeError > 1E-14) {
        teststatus = false;
    }
    
    msm4g_simulation_delete(simulation);
    return teststatus;
}
