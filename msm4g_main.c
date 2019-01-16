/** @file msm4g_main.c
 * @brief The main entry point of the MSM4G application.
 *
 * @mainpage
 *
 * MSM4G program is aimed to use the Multi-Level Summation Method (MSM) for (4)
 * Gravitational (G) N-Body problems @anchor logonote *.
 *
 * @section datastructures Data Structures
 *
 * - Particle is used to represent a celestial object.
 * - LinkedList is a generic two-way linked listed implementation.
 * - Bin is a container to hold the particles in geometric hashing.
 * - SimulationBox is a rectangular box enclosing all of the particles.
 * - SimulationParameters contains relative cutoff, order of B-splines, etc.
 *
 * @section shortrange Short-Range Calculations
 * 
 * The simulation domain is divided into cubic boxes called Bin whose edges are
 * set to the cut-off parameter. Depending on the particle distribution and the 
 * arrangement of the simulation box, there might be considerable amount of
 * Bin's containing no particles. In order to save some space, the Bin's are
 * stored as a linked list. Only the Bin's containing at least one particle are
 * stored in the linked list.
 * 
 * The neighbors of a specific Bin are also stored in a linked list. Again,
 * a Bin is added as a neighbor only if it contains at least a particle.
 *
 * In order to exploit the Newton's third law during the short-range
 * calculations, each Bin is enumerated by a unique integer corresponding
 * to its index in the three dimensional simulation box. An exotic choise is
 * Cantor's generalized pairing function (See. msm4g_math_cantor()). Whenever
 * a new Bin is allocated, its cantor index is calculated and stored in the
 * Bin structure. Since the Bin structure already has the index information,
 * there is no need to invert the cantor index.
 *
 * @section tests Test Functions
 * 
 * A user who wants to use MSM4G in his/her application or
 * tries to understand the details of the package may prefer
 * to look at the test functions in msm4g_tests.c as they are
 * self-contained and test various features individually.
 *
 * @ref logonote "*" The artwork at top of the page is the picture of the
 * Andromeda Galaxy which is taken using amateur refractor telescope.
 * Courtesy of Kanwar Singh. The link to the higher resolutions is:
 * https://commons.wikimedia.org/wiki/File:ANDROMEDA_GALAXY.jpg
 */

#include "msm4g_lib.h"

int main(int argc,char *argv[]) {
  if ( argc != 12 ) {
    fprintf(stderr, "Usage: %s datafile Ax Ay Az abar mu nu Mx My Mz L\n",
        argv[0]);
    fprintf(stderr, "Set Mx=My=Mz=L=0 to switch automatic selection\n");
    fprintf(stderr,
        "Example: ./msm4g data/changaN512.ini 1 1 1 4 2 4 0 0 0 0\n");
    return 1;
  }

  Boolean periodic = true;


  SimulationBox *box = msm4g_box_new();
  box->x = 0;
  box->y = 0;
  box->z = 0;
  box->wx = atof(argv[2]);
  box->wy = atof(argv[3]);
  box->wz = atof(argv[4]);
  double abar = atof(argv[5]);
  int mu = atoi(argv[6]);
  int nu = atoi(argv[7]);
  int Mx = atoi(argv[8]);
  int My = atoi(argv[9]);
  int Mz = atoi(argv[10]);
  int L =  atoi(argv[11]);
  Simulation *simulation = msm4g_simulation_new(argv[1], box, periodic, nu,
      abar,mu,L,Mx,My,Mz);

  msm4g_simulation_run(simulation);

  int N = simulation->parameters->N ;

  FILE *fp=fopen("msm.acc","w");
  fprintf(fp,"%d\n",N);
  for (int i = 0 ; i < N ; i++) {
    double totalx = simulation->particles[i].acc_total[0];
    fprintf(fp,"%25.16e\n",totalx);
  }
  for (int i = 0 ; i < N ; i++) {
    double totaly = simulation->particles[i].acc_total[1];
    fprintf(fp,"%25.16e\n",totaly);
  }
  for (int i = 0 ; i < N ; i++) {
    double totalz = simulation->particles[i].acc_total[2];
    fprintf(fp,"%25.16e\n",totalz);
  }
  fclose(fp);

#ifdef DEBUG
  fp=fopen("msm.acc_noninterpolated","w");
  fprintf(fp,"%d\n",N);
  for (int i = 0; i < N ;i++) {
    Particle *p = &(simulation->particles[i]);
    int j = 0;
    fprintf(fp,"%25.16e %25.16e %25.16e\n",p->acc_short_true[j]+p->acc_long_direct_true[j],
           p->acc_long_fourier_true[j],
           p->acc_short_true[j]+p->acc_long_direct_true[j]+p->acc_long_fourier_true[j]);
  }
  for (int i = 0; i < N ;i++) {
    Particle *p = &(simulation->particles[i]);
    int j = 1;
    fprintf(fp,"%25.16e %25.16e %25.16e\n",p->acc_short_true[j]+p->acc_long_direct_true[j],
           p->acc_long_fourier_true[j],
           p->acc_short_true[j]+p->acc_long_direct_true[j]+p->acc_long_fourier_true[j]);
  }
  for (int i = 0; i < N ;i++) {
    Particle *p = &(simulation->particles[i]);
    int j = 2;
    fprintf(fp,"%25.16e %25.16e %25.16e\n",p->acc_short_true[j]+p->acc_long_direct_true[j],
           p->acc_long_fourier_true[j],
           p->acc_short_true[j]+p->acc_long_direct_true[j]+p->acc_long_fourier_true[j]);
  }
  fclose(fp);
#endif
  
  fp = fopen("msm.pot","w");
  fprintf(fp,"%25.16e\n",simulation->output->potentialEnergyTotal);
  fclose(fp);
  msm4g_simulation_save(simulation, stdout);
  msm4g_simulation_delete(simulation);

  return 0;
}
