/** @file msm4g_tests.c
 * @brief The test functions to check the functionality of the MSM4G package.
 *
 * This is probably the best place one can learn how to use the package.
 * These test functions are used to check the outcomes of the functions
 * by comparing with expected results.
 */
#include "msm4g_lib.h"
#include "msm4g_tests.h"

char    testnames[100][64];
Boolean teststatus[100];
int     testcount  = 0 ;

int main() {
  msm4g_unit_test_all();
}
void msm4g_test_summary() {
  int failed = 0;
  for (int i = 0 ; i < testcount ; i++) {
    printf("%02d Testing %-64s : ",i, testnames[i]);
    if (teststatus[i]) {
      printf("passed\n");
    } else {
      printf("failed\n");
      failed++;
    }
  }
  printf("%2d total, %2d passed, %d failed\n",
      testcount,testcount-failed,failed);
}

void msm4g_test_assert(const char *name,int status) {
  strcpy(testnames[testcount],name);
  teststatus[testcount] = status;
  testcount++;
  if (testcount >= 100) {
    fprintf(stderr,
        "Too many test cases, please increase the size of test slot\n");
  }
}

void msm4g_unit_test_all() {
  int numberoftests;
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
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_16);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_17);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_18);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_19);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_20);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_21);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_22);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_23);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_24);
  msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_25);

  numberoftests = msm4g_linkedlist_size(list);

  for (int i = 0; i < numberoftests ; i++) {
    testFunction = (testFunctionType *)msm4g_linkedlist_get(list, i);
    testFunction();
  }

  msm4g_test_summary();

  msm4g_linkedlist_destroy(list);

}

Boolean msm4g_unit_test_1() {
  Boolean status = true;
  LinkedList *list;
  Particle x[10];
  int i;

  list = msm4g_linkedlist_new();

  for (i=0;i<10;i++) {
    msm4g_particle_reset(x+i);
    msm4g_linkedlist_add(list,&x[i]);
  }
  if (msm4g_linkedlist_size(list) != 10) {
    msm4g_linkedlist_destroy(list);
    status = false;
  }
  msm4g_linkedlist_destroy(list);
  msm4g_test_assert("Counting the members of a two-way linked list",
      status == true);
  return status;
}

Boolean msm4g_unit_test_2() {
  LinkedList *list;
  LinkedListElement *curr;
  Particle x[10];
  Particle *y;
  Boolean status ;
  int i;

  status = true;
  list = msm4g_linkedlist_new();

  for (i=0;i<10;i++) {
    x[i].m = i;
    msm4g_linkedlist_add(list,&x[i]);
  }

  /**< Iterating from tail to head */
  i=9;
  curr=list->tail;
  while (curr != NULL) {
    y =  (Particle *)curr->data;
    if (fabs(i - (y->m)) > DBL_EPSILON) {
      status = false;
      break;
    }
    curr = curr->prev;
    i--;
  }
  msm4g_linkedlist_destroy(list);
  msm4g_test_assert("Traversing a linked list from tail to head",
      status == true);
  return status;
}

Boolean msm4g_unit_test_3() {
  LinkedList *list;
  LinkedListElement *curr;
  Particle x[10];
  Particle *y;
  Boolean status ;
  int i;

  status = true;
  list = msm4g_linkedlist_new();

  for (i=0;i<10;i++) {
    x[i].m = i;
    msm4g_linkedlist_add(list,&x[i]);
  }

  i=0;
  curr=list->head;
  while (curr != NULL) {
    y =  (Particle *)curr->data;
    if (fabs(i - (y->m)) > DBL_EPSILON) {
      status = false;
      break;
    }
    curr = curr->next;
    i++;
  }
  msm4g_linkedlist_destroy(list);
  msm4g_test_assert("Traversing a linked list from head to tail",
      status == true);
  return status;
}


Boolean msm4g_unit_test_4() {
  Boolean status = true;
  const int DIM = 3;
  const int N = 3; /* Number of particles in eight.ini */
  int dummy;
  Particle particles[3];
  Particle *particlesInFile;
  int i,j;

  for (i=0;i<N;i++) {
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


  for (i=0;i<N;i++) {
    if (fabs(particles[i].m - particlesInFile[i].m)>DBL_EPSILON ) {
      status = false;
      break;
    }
    for (j=0;j<DIM;j++) {
      if (fabs(particles[i].r.value[j]-particlesInFile[i].r.value[j])
          > DBL_EPSILON) {
        status = false;
        break;
      }
      if (fabs(particles[i].v.value[j]-particlesInFile[i].v.value[j])
          > DBL_EPSILON) {
        status = false;
        break;
      }
    }
  }

  /* Deallocate the particles */
  free(particlesInFile);
  msm4g_test_assert("Reading particles from text file", status == true);
  return status;
}

Boolean msm4g_unit_test_5() {
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
  for (i=0;i<n;i++) {
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
  p=2; /* degree of base polynomials.
  if p=1 there is no need for padding around the boundary */
  msm4g_box_update(box, list, 0.5,h,p);

  msm4g_d3vector_daxpy(&locationError, &locationExpected,-1.0,&(box->location));
  msm4g_d3vector_daxpy(&widthError,    &widthExpected,   -1.0,&(box->width));

  if (msm4g_d3vector_norm(&locationError) > DBL_EPSILON ||
      msm4g_d3vector_norm(&widthError)    > DBL_EPSILON ) {
    status = false;
  }

  msm4g_linkedlist_destroy(list);
  msm4g_particle_destroyarray(particles, n);
  msm4g_box_destroy(box);
  msm4g_test_assert("Simulation box create and update", status == true);
  return status;
}

Boolean msm4g_unit_test_6() {
  Boolean status = true;
  Boolean periodic = false;
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
  binlist=msm4g_bin_generate(&box,particlelist,6,binwidth,periodic);

  bin = (Bin *)msm4g_linkedlist_get(binlist, 0);
  if (bin->cantorindex != 47) status = false;
  /* It should has only one neighbor */
  if (msm4g_linkedlist_size(bin->neighbors) != 1) status = false;
  /* Its neighbor's cantor index should be 17 */
  if (((Bin *)msm4g_linkedlist_get(bin->neighbors,0))->cantorindex != 73)
    status = false;
  /* It should contain two particles */
  if (msm4g_linkedlist_size(bin->particles) != 2) status = false;

  bin = (Bin *)msm4g_linkedlist_get(binlist, 1);
  /* This bin should has a neighbor */
  if (msm4g_linkedlist_size(bin->neighbors) != 0) status = false;
  /* Its cantor index should be 97 */
  if (bin->cantorindex != 97) status = false;
  /* It should contain only one particles */
  if (msm4g_linkedlist_size(bin->particles) != 1) status = false;

  bin = (Bin *)msm4g_linkedlist_get(binlist, 2);
  /* It should has only one neighbor */
  if (msm4g_linkedlist_size(bin->neighbors) != 1) status = false;
  /* Its neighbor's cantor index should be 8 */
  if (((Bin *)msm4g_linkedlist_get(bin->neighbors,0))->cantorindex != 47)
    status = false;
  /* It should contain three particles */
  if (msm4g_linkedlist_size(bin->particles) != 3) status = false;

  msm4g_bin_destroy(binlist);
  free(particlelist);
  msm4g_test_assert("Binning the simulation box for short-range calculation",
      status == true);

  return status;
}

Boolean msm4g_unit_test_7() {
  Boolean status = true;
  int i,j;
  int inp[6] = {0, 1, 2, 3,  4, 5  };
  int out[6] = {1, 1, 2, 6, 24, 120};
  int x[4][4];
  int xexpected[4][4] = {{0,1,3,6},{2,4,7,11},{5,8,12,17},{9,13,18,24}};
  int vector2d[2];

  /* Checking factorial function */
  for (i=0; i<6; i++)
    if (msm4g_math_factorial(inp[i]) != out[i]) {
      status = false;
      break;
    }
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      vector2d[0] = i; vector2d[1]=j;
      x[i][j] = msm4g_math_cantor(vector2d,2);
      if (x[i][j] != xexpected[i][j]) {
        status = false;
        break;
      }
    }
  }
  msm4g_test_assert("Cantor pairing function", status == true);
  return status;
}

Boolean msm4g_unit_test_8() {
  Boolean status = true;
  int i;
  int n = 3;

  for (i=0; i<n; i++) {
    /* \gamma(1.0)  = 1.0 for all gamma function. */
    if (fabs(msm4g_smoothing_gama(1,i+1)-1.0) > DBL_EPSILON) {
      status = false;
      break;
    }
    /* \gamma'(1.0) = -1.0 for all gamma function. */
    if (fabs(msm4g_smoothing_gamaprime(1,i+1)+1.0) > DBL_EPSILON) {
      status = false;
      break;
    }
    /* \gamma'(0.0) = 0.0 for all gamma function. */
    if (fabs(msm4g_smoothing_gamaprime(0,i+1)-0.0) > DBL_EPSILON) {
      status = false;
      break;
    }

  }
  /* \gamma(0.0) = 3/2 for C1 gamma */
  if (fabs(msm4g_smoothing_gama(0,2)-3.0/2.0)   > DBL_EPSILON) {
    status = false;
  }
  /* \gamma(0.0) = 15/8 for C2 gamma */
  if (fabs(msm4g_smoothing_gama(0,3)-15.0/8.0)  > DBL_EPSILON) {
    status = false;
  }
  /* \gamma(0.0) = 35/16 for C3 gamma */
  if (fabs(msm4g_smoothing_gama(0,4)-35.0/16.0) > DBL_EPSILON) {
    status = false;
  }

  msm4g_test_assert("Smoothing functions aka even-powered softeners",
      status == true);
  return status;
}

Boolean msm4g_unit_test_9() {
  Boolean teststatus = true;
  Boolean periodic = false;
  int order = 2 ;
  double abar = 1;
  int mu = 2;
  char *datafile = "data/bintest.ini";
  double boxlocation = 0.0;
  double boxwidth = 30.0;
  SimulationBox *box = msm4g_box_newCube(boxlocation,boxwidth);
  Simulation *simulation = msm4g_simulation_new(datafile,box,periodic,order,
      abar,mu,0,0,0,0);
  simulation->parameters->a = 10;
  SimulationParameters *sp = simulation->parameters;
  Particle *particles = simulation->particles;
  LinkedList *binlist = msm4g_bin_generate(box,particles,
      simulation->parameters->N,sp->a,periodic);
  msm4g_force_short(binlist, sp->a, simulation);
  /* Calculate short-range potential energy */
  double energy = 0.0;
  for (int i = 0 ; i < simulation->parameters->N ; i++) {
    energy += particles[i].potential_short_real * particles[i].m ;
  }
  simulation->output->potentialEnergyShortRange = energy;

  msm4g_bin_destroy(binlist);

  double potential = simulation->output->potentialEnergyShortRange ;
  double potentialExpected =  1127.0/6000 + (5000*sqrt(2.0)-1776.0)/4000.0 +
      (2000*sqrt(82.0)-17876.0)/164000.0;

  if (fabs(potentialExpected-potential)/potential  > 1E-14) {
    teststatus = false;
  }

  msm4g_simulation_delete(simulation);
  msm4g_test_assert("Short-range potential energy for a 6-particle case",
      teststatus == true);

  return teststatus;
}

Boolean msm4g_unit_test_10() {
  Boolean status = true;
  AbstractGrid *grid;

  double h = 1.0;
  int nx = 3;
  int ny = 3;
  int nz = 3;

  grid = msm4g_grid_dense_new(nx,ny,nz,h,h,h);
  /* Check if it could allocated the object */
  if (grid == NULL) status = false;

  grid->reset(grid,0.0);

  grid->setElement(grid,0,1,2,10.0);

  msm4g_grid_destroy(&grid);
  /* The grid should point to NULL after destruction */
  if (grid != NULL) status = false;

  msm4g_test_assert("Dense grid implementation", status == true);

  return status;
}

Boolean msm4g_unit_test_11() {
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
  Simulation *simulation = msm4g_simulation_new(datafile,box,periodic,4,
      abar,mu,0,0,0,0);
  simulation->parameters->Mx = 1;
  simulation->parameters->My = 1;
  simulation->parameters->Mz = 1;
  simulation->parameters->h  = h ;

  msm4g_anterpolation(simulation);

  grid = simulation->gridmass[0];
  for (int i = 0 ; i < grid->nx ; i++)
    for (int j = 0 ; j < grid->ny ; j++)
      for (int k = 0 ; k < grid->nz ; k++)
        if (fabs(grid->getElement(grid,i,j,k)-phix[i]*phix[j]*phix[k])
            > FLT_EPSILON) {
          teststatus = false;
          break;
        }


  msm4g_simulation_delete(simulation);
  msm4g_test_assert("Anterpolation for a special single-particle case",
      teststatus == true);
  return teststatus;
}

Boolean msm4g_unit_test_12() {
  Boolean teststatus = true;
  /* Expected values are from Table I of doi:10.1063/1.4943868 */
  char *omegaprimeCubicExpected[13] = {"3.464","-1.732", "0.679","-0.240",
      "0.080","-0.026","0.008","-0.002","0.001","-0.000",
      "0.000","-0.000","0.000"};
  char *omegaprimeQuinticExpected[13] = {"12.379","-9.377", "5.809","-3.266",
      "1.735","-0.889","0.444","-0.217","0.105","-0.050",
      "0.024","-0.011","0.005"};
  double *omegaprimeCubicCalculated = msm4g_util_omegaprime(20, 4);
  double *omegaprimeQuinticCalculated = msm4g_util_omegaprime(20, 6);
  for (int i = 0 ; i < 13 ; i++) {
    char *expectedCubic   = omegaprimeCubicExpected[i];
    char *expectedQuintic = omegaprimeQuinticExpected[i];
    char calculated[7]  ;
    sprintf(&(calculated[0]),"%5.3f",omegaprimeCubicCalculated[i]);
    if (strcmp(expectedCubic, calculated) != 0) {
      teststatus = false;
      break;
    }
    sprintf(&(calculated[0]),"%5.3f",omegaprimeQuinticCalculated[i]);
    if (strcmp(expectedQuintic, calculated) != 0) {
      teststatus = false;
      break;
    }
  }
  free(omegaprimeCubicCalculated);
  free(omegaprimeQuinticCalculated);
  msm4g_test_assert("Quasi-interpolation coefficients (omegaprime)",
      teststatus == true);
  return teststatus;
}

Boolean msm4g_unit_test_13() {
  Boolean teststatus = true;
  const int ntest = 9;
  const int cubic = 4;
  const int quintic = 6;
  const double mathematicaCubic[9]        =
  {0 ,1,   8,   23,   32,   23,   8,  1,0}; /* divided by    48 */
  const double mathematicaCubicPrime[9]   =
  {0 ,1,   4,    5,    0,   -5,  -4, -1,0}; /* divided by     8 */
  const double mathematicaQuintic[9]      =
  {0,81,2528,13438,22528,13438,2528, 81,0}; /* divided by 40960 */
  const double mathematicaQuinticPrime[9] =
  {0,27, 400,  942,    0, -942,-400,-27,0}; /* divided by  2048 */

  for (int i = 0 ; i < ntest ; i++) {
    double expectedCubic        = mathematicaCubic[i]       /   48.0 ;
    double expectedCubicPrime   = mathematicaCubicPrime[i]  /    8.0 ;
    double expectedQuintic      = mathematicaQuintic[i]     /40960.0 ;
    double expectedQuinticPrime = mathematicaQuinticPrime[i]/ 2048.0 ;
    double calculatedCubic        = msm4g_bases_bspline     (cubic  ,
        cubic  *i/(double)(ntest-1));
    double calculatedCubicPrime   = msm4g_bases_bsplineprime(cubic  ,
        cubic  *i/(double)(ntest-1));
    double calculatedQuintic      = msm4g_bases_bspline     (quintic,
        quintic*i/(double)(ntest-1));
    double calculatedQuinticPrime = msm4g_bases_bsplineprime(quintic,
        quintic*i/(double)(ntest-1));
    if (fabs(expectedCubic       -calculatedCubic)        > DBL_EPSILON) {
      teststatus = true;
      break;
    }
    if (fabs(expectedCubicPrime  -calculatedCubicPrime)   > DBL_EPSILON) {
      teststatus = true;
      break;
    }
    if (fabs(expectedQuintic     -calculatedQuintic)      > DBL_EPSILON) {
      teststatus = true;
      break;
    }
    if (fabs(expectedQuinticPrime-calculatedQuinticPrime) > DBL_EPSILON) {
      teststatus = true;
      break;
    }
  }
  msm4g_test_assert("Implementation of B-Spline and its derivatives",
      teststatus == true);
  return teststatus;
}

Boolean msm4g_unit_test_14() {
  Boolean teststatus = true;
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
      if (fabs(gama     -expectedgama[nu-2][i])      > DBL_EPSILON) {
        teststatus = false;
        break;
      }
      if (fabs(gamaprime-expectedgamaprime[nu-2][i]) > DBL_EPSILON) {
        teststatus = false;
        break;
      }
    }
  }
  msm4g_test_assert("Even-powered softener", teststatus == true);
  return teststatus;
}

Boolean msm4g_unit_test_15() {
  char message[100];
  Boolean teststatus = true;
  Boolean periodic = true;
  int order = 4;
  int mu = 2;
  double abar = 4.0;
  double expected = 6.9843177228211383e-09;
  double expectedacc[8][3] = {
    { 1.7127138343151520e-05,-7.9922393447316331e-06, 7.0828748251281019e-06},
    {-3.1966280586136999e-05,-1.0376534820290217e-04,-6.0893167284709774e-05},
    {-4.3875970537235183e-05,-2.7824718658855566e-05, 3.0876280003449423e-05},
    { 7.3727156274532452e-05, 1.1860770584138265e-04,-1.1113166484930614e-05},
    { 7.8059074802170574e-06,-1.4019342701374309e-07,-5.7782346576733895e-05},
    { 5.9620212096497241e-05,-2.3746627592581602e-05, 3.1717687081769986e-05},
    { 1.4998005198020302e-05, 2.9271945320370181e-05, 6.7963357688293238e-06},
    {-9.7436168269046449e-05, 1.5589476064332005e-05, 5.3315502667197417e-05}};

  SimulationBox *unitCube = msm4g_box_newCube(0,1);

  Simulation *simulation = msm4g_simulation_new("data/changaN8.ini",
      unitCube,periodic,order,abar,mu,0,0,0,0);

  SimulationParameters *sp = simulation->parameters;
  SimulationBox *box = simulation->box;
  Particle *particles = simulation->particles;

  LinkedList *binlist = msm4g_bin_generate(box,particles,
      simulation->parameters->N,sp->a,periodic);

  msm4g_force_short(binlist, sp->a, simulation);

  /* Calculate short-range potential energy */
  { double energy = 0.0;
    for (int i = 0 ; i < simulation->parameters->N ; i++) {
      energy += particles[i].potential_short_real * particles[i].m ;
    }
    simulation->output->potentialEnergyShortRange = energy;
  }

  msm4g_bin_destroy(binlist);

  double relativeError = fabs(simulation->output->potentialEnergyShortRange-
      expected)/fabs(expected);
  sprintf(message,"Short-range potential energy for ChaNGaN8: %8.2e",
      relativeError);
  msm4g_test_assert(message,relativeError < 1E-14);

  if (relativeError > 1E-14) {
    teststatus = false;
  }

  for (int i = 0 ; i < simulation->parameters->N ; i++) {
    double relerr = msm4g_util_diffnorm(expectedacc[i],
        simulation->particles[i].acc_short, 3) /
        msm4g_util_norm(expectedacc[i],3);
    if (relerr > 1E-14) {
      teststatus = false;
      break;
    }
  }
  msm4g_test_assert("Short-range force for ChaNGa N=8",  teststatus == true );

  msm4g_simulation_delete(simulation);
  return teststatus;
}

Boolean msm4g_unit_test_16() {
  Boolean teststatus = true;
  Boolean periodic = true;
  int order = 4;
  int mu = 2;
  double abar = 4.0;

  double expectedGridMass[8] = {1.0634562200737233e-05,
      1.0640213996318175e-05,
      9.5181244646372011e-06,
      9.8543177848436609e-06,
      1.0582389591755571e-05,
      1.0377630551436881e-05,
      9.0976320916671066e-06,
      9.2398733186041778e-06};
  SimulationBox *unitCube = msm4g_box_newCube(0,1);
  Simulation *simulation = msm4g_simulation_new("data/changaN8.ini",
      unitCube,periodic,order,abar,mu,0,0,0,0);

  msm4g_anterpolation(simulation);

  AbstractGrid *grid = simulation->gridmass[0] ;
  int counter = 0;
  for (int mx = 0 ; mx < simulation->parameters->Mx ; mx++) {
    for (int my = 0 ; my < simulation->parameters->My ; my++) {
      for (int mz = 0 ; mz < simulation->parameters->Mz ; mz++) {
        double expected = expectedGridMass[counter++] ;
        double calculated = grid->getElement(grid,mx,my,mz);
        double relerr = fabs(expected-calculated)/fabs(expected);
        if (relerr > 1E-14) {
          teststatus = false;
          break;
        }
      }
    }
  }
  msm4g_test_assert("Anterpolation for ChaNGa N=8",  teststatus == true );
  msm4g_simulation_delete(simulation);
  return teststatus;
}

Boolean msm4g_unit_test_17() {
  Boolean teststatus = true;
  double h0 = 1./4 ;
  for (int i = 0 ; i < 80 ; i++) {
    double aL = 2 + i/10.0 ;
    double beta = msm4g_util_choose_beta(aL,TOL_DIRECT,h0);
    double res = fabs(erfc(beta*aL) - TOL_DIRECT * aL / h0) ;
    if (res > 1E-13) {
      teststatus = false;
      break;
    }
    double kmax = msm4g_util_choose_kmax(beta, TOL_FOURIER, h0);
    res = fabs(erfc(MYPI*kmax/beta) - sqrt(MYPI)*TOL_FOURIER/(2*beta*h0));
    if (res > 1E-13) {
      teststatus = false;
      break;
    }
  }
  msm4g_test_assert("Choosing optimal Ewald's splitting parameter",
      teststatus == true );
  return teststatus;
}

Boolean msm4g_unit_test_18() {
  Boolean teststatus = true;
  Simulation *simulation;
  char message[100];
  SimulationBox *box = msm4g_box_newCube(0, 1);
  double abar = 4; int mu = 6, nu=4;

  simulation = msm4g_simulation_new("data/CsClN2.ini", box, true, nu,
      abar, mu,0,0,0,0);
  msm4g_simulation_run(simulation);
  double calculated = simulation->output->potentialEnergyTotal;
  double expected = -2.035361508229;
  double relerr=fabs(calculated-expected)/fabs(expected);
  sprintf(message,"Rel. err. in pot. energy of CsCl N=2 is  %9.2E<1E-3",relerr);
  msm4g_test_assert(message,relerr < 1E-3);
  msm4g_simulation_delete(simulation);
  return teststatus ;
}

Boolean msm4g_unit_test_19() {
  Boolean teststatus = true;
  Simulation *simulation;
  char message[100];
  SimulationBox *box = msm4g_box_newCube(0, 2);
  double abar = 4; int mu = 6, nu=4;
  
  simulation = msm4g_simulation_new("data/CsClN16.ini", box, true, nu,
                                    abar, mu,0,0,0,0);
  msm4g_simulation_run(simulation);
  double calculated = simulation->output->potentialEnergyTotal;
  double expected = -2.035361508229 * 8;
  double relerr=fabs(calculated-expected)/fabs(expected);
  sprintf(message,"Rel. err. in pot. energy of CsCl N=16 is  %9.2E<1E-3",relerr);
  msm4g_test_assert(message,relerr < 1E-3);
  msm4g_simulation_delete(simulation);
  return teststatus ;
}

Boolean msm4g_unit_test_20() {
  Boolean teststatus = true;
  double expectedValues[2][5] = {{ 6/ 8. ,  4/ 8.0, 1/ 8.0, 0     , 0},
      {20/32.0, 15/32.0, 6/32.0, 1/32.0, 0}};
  for (int v = 4 ; v <= 6; v += 2) {
    for (int k = 0 ; k < 5 ;k++) {
      double expected = expectedValues[v/2-2][k] ;
      double calculated = msm4g_util_jn(v, k);
      double relerr = fabs(calculated-expected)/fabs(expected);
      if (relerr > 1E-14) {
        teststatus = false ;
        break;
      }
    }
  }
  msm4g_test_assert("JN utility function", teststatus == true);
  return teststatus;
}

Boolean msm4g_unit_test_21() {
  Boolean teststatus = true;
  Boolean periodic = true;
  int mu = 2;
  int nu = 4;
  double abar = 4;
  char message[100] ;
  SimulationBox *box = msm4g_box_newCube(0, 2);
  Simulation *simulation = msm4g_simulation_new("data/NaClN8.ini", box,
      periodic, nu, abar,mu,0,0,0,0);
  int N = simulation->parameters->N ;
  int L = simulation->parameters->L ;
  msm4g_simulation_run(simulation);
  double qsum = 0.0;
  for (int i = 0 ; i < N ; i++)
    qsum += simulation->particles[i].m ;
  for (int l = 0 ; l <= L ; l++ ) {
    double q1sum = msm4g_grid_dense_sum(simulation->gridmass[l]);
    sprintf(message,"NaClN8: sum of grid masses at l=%d stays same",l);
    msm4g_test_assert(message,  fabs(qsum-q1sum) < 1E-14);
  }
  double energy = simulation->output->potentialEnergyTotal ;
  double energyExpected = - 1.747564594633182 * N * 0.5 ;
  double relerr = fabs(energy-energyExpected)/fabs(energyExpected);
  sprintf(message,"Rel. err. in pot. energy of NaCl N=8 is %9.2e<1E-3",relerr);
  msm4g_test_assert(message,relerr < 1E-3);
  msm4g_simulation_delete(simulation);
  return teststatus ;
}

Boolean msm4g_unit_test_22() {
  Boolean teststatus = true;
  Boolean periodic = true;
  int mu = 2;
  int nu = 4;
  double abar = 4;
  char message[100] ;
  SimulationBox *box = msm4g_box_newCube(0, 4);
  Simulation *simulation = msm4g_simulation_new("data/NaClN64.ini", box,
      periodic, nu, abar,mu,0,0,0,0);
  int N = simulation->parameters->N ;
  int L = simulation->parameters->L ;
  msm4g_simulation_run(simulation);
  double qsum = 0.0;
  for (int i = 0 ; i < N ; i++)
    qsum += simulation->particles[i].m ;
  for (int l = 0 ; l <= L ; l++ ) {
    double q1sum = msm4g_grid_dense_sum(simulation->gridmass[l]);
    sprintf(message,"NaClN64: sum of grid masses at l=%d stays same",l);
    msm4g_test_assert(message,  fabs(qsum-q1sum) < 1E-14);
  }
  double energy = simulation->output->potentialEnergyTotal ;
  double energyExpected = - 1.747564594633182 * N * 0.5 ;
  double relerr = fabs(energy-energyExpected)/fabs(energyExpected);
  sprintf(message,"Rel. err. in pot. energy of NaCl N=64 is %9.2e<1E-3",relerr);
  msm4g_simulation_delete(simulation);
  return teststatus ;
}

Boolean msm4g_unit_test_23() {
  Boolean teststatus = true;
  Boolean periodic = true;
  int mu = 2;
  int nu = 4;
  double abar = 4;
  char message[100];
  SimulationBox *box = msm4g_box_newCube(0, 8);
  Simulation *simulation = msm4g_simulation_new("data/NaClN512.ini", box,
      periodic, nu, abar,mu,0,0,0,0);
  int N = simulation->parameters->N ;
  int L = simulation->parameters->L ;
  msm4g_simulation_run(simulation);
  double qsum = 0.0;
  for (int i = 0 ; i < N ; i++)
    qsum += simulation->particles[i].m ;
  for (int l = 0 ; l <= L ; l++ ) {
    double q1sum = msm4g_grid_dense_sum(simulation->gridmass[l]);
    sprintf(message,"NaClN512: sum of grid masses at l=%d stays same",l);
    msm4g_test_assert(message,  fabs(qsum-q1sum) < 1E-14);
  }

  double energy = simulation->output->potentialEnergyTotal ;
  double energyExpected = - 1.747564594633182 * N * 0.5 ;
  double relerr = fabs(energy-energyExpected)/fabs(energyExpected);
  sprintf(message,"Rel. err. in pot. energy of NaCl N=512 is %9.2e<1E-3"
      ,relerr);
  msm4g_test_assert(message,relerr < 1E-3);
  msm4g_simulation_delete(simulation);
  return teststatus ;
}

Boolean msm4g_unit_test_24() {
  Boolean teststatus = true;
  Simulation *simulation;
  char message[100];
  SimulationBox *box = msm4g_box_newCube(0, 4);
  double abar = 4; int mu = 6, nu=4;
  
  simulation = msm4g_simulation_new("data/CsClN128.ini", box, true, nu,
                                    abar, mu,0,0,0,0);
  msm4g_simulation_run(simulation);
  double calculated = simulation->output->potentialEnergyTotal;
  double expected = -2.035361508229 * 8 * 8;
  double relerr=fabs(calculated-expected)/fabs(expected);
  sprintf(message,"Rel. err. in pot. energy of CsCl N=128 is  %9.2E<1E-3",relerr);
  msm4g_test_assert(message,relerr < 1E-3);
  msm4g_simulation_delete(simulation);
  return teststatus ;
}

Boolean msm4g_unit_test_25() {
  Boolean teststatus = true;
  Simulation *simulation;
  char message[100];
  SimulationBox *box = msm4g_box_newCube(0, 1);
  double abar = 4; int mu2 = 2, mu6 = 6, nu=4;
  
  simulation = msm4g_simulation_new("data/CsClN2.ini", box, true, nu,
                                    abar, mu2,0,0,0,0);
  msm4g_simulation_run(simulation);
  double calculated = simulation->output->potentialEnergyTotal;
  double expected = -2.035361508229;
  double relerr=fabs(calculated-expected)/fabs(expected);
  msm4g_simulation_delete(simulation);
  
  box = msm4g_box_newCube(0, 1);
  simulation = msm4g_simulation_new("data/CsClN2.ini", box, true, nu,
                                    abar, mu6,0,0,0,0);
  msm4g_simulation_run(simulation);
  calculated = simulation->output->potentialEnergyTotal;
  expected = -2.035361508229;
  double relerr2=fabs(calculated-expected)/fabs(expected);
  msm4g_simulation_delete(simulation);
  
  sprintf(message,"Increasing mu=2 to 6 decrease error CsCl N=2 %9.2E<%9.2E",relerr2,relerr);
  msm4g_test_assert(message,relerr2 < relerr);
  return teststatus ;
}
