/** @file msm4g_lib.c
 * @brief The definitions of all core functions of MSM4G package.
 */
#include "msm4g_lib.h"
#include "msm4g_bases.h"

void msm4g_tic() {
  msm4g_tictocmanager(1);
}

double msm4g_toc() {
  return msm4g_tictocmanager(0);
}

double msm4g_tictocmanager(int push) {
  double elapsed_seconds = 0.0;
  static clock_t stack_data[100] ;
  static int stack_lastindex = 0 ;
  if (push) {
    stack_data[stack_lastindex] = clock();
    stack_lastindex = stack_lastindex + 1;
  } else {
    clock_t now = clock();
    stack_lastindex = stack_lastindex - 1;
    clock_t previous = stack_data[stack_lastindex];
    elapsed_seconds = (double)(now-previous)/CLOCKS_PER_SEC;
  }
  return elapsed_seconds;
}

double msm4g_potential_energy(Simulation *simulation) {
  double beta = simulation->parameters->beta ;
  double a = simulation->parameters->a ;
  double Ax = simulation->box->wx;
  double Ay = simulation->box->wy;
  double Az = simulation->box->wz;
  double detA = Ax * Ay * Az ;
  int L = simulation->parameters->L ;
  int N = simulation->parameters->N;
  int nu = simulation->parameters->nu;
  Particle *particles = simulation->particles;
  SimulationOutput *so = simulation->output;

  double ushort_real = 0.0;
  for (int i = 0 ; i < N ; i++) {
    ushort_real += particles[i].potential_short_real * particles[i].m ;
  }

#ifdef DEBUG
  double ushort_real_true = 0.0;
  for (int i = 0 ; i < N ; i++) {
    ushort_real_true += particles[i].potential_short_real_true * particles[i].m ;
  }
  so->ushort_real_true = ushort_real_true;
#endif

  double csr = MYPI / (beta * beta * detA);
  double qsum = 0.0;
  for (int i = 0; i < N; i++)
    qsum += particles[i].m;
  double ushort_csr = - 0.5 * qsum * qsum * csr;

  double ushort_self = 0.0;
  int pxmin = - a / Ax;
  int pxmax =   a / Ax;
  int pymin = - a / Ay;
  int pymax =   a / Ay;
  int pzmin = - a / Az;
  int pzmax =   a / Az;
  double psum = 0.0;
  for (int px = pxmin; px <= pxmax; px++) {
    for (int py = pymin; py <= pymax; py++) {
      for (int pz = pzmin; pz <= pzmax; pz++) {
        if (px == 0 && py == 0 && pz == 0) continue;
        double rlen2 = Ax * px * Ax * px
            + Ay * py * Ay * py
            + Az * pz * Az * pz;
        double rlen = sqrt(rlen2);
        psum += msm4g_kernel(0,L,rlen,a,beta,nu);
      }
    }
  }
  for (int i = 0; i < N; i++) {
    ushort_self += 0.5 * particles[i].m * particles[i].m * psum ;
  }


  double q2sum = 0.0;
  for (int i = 0; i < N; i++)
    q2sum += particles[i].m * particles[i].m;
  double ulong_self = - 0.5 * q2sum * (1.0 / a) * msm4g_smoothing_gama(0, nu);

  /* Calculate long-range potential energy */
  AbstractGrid *e1 = simulation->gridpotential[0];
  AbstractGrid *q1 = simulation->gridmass[0];
  double ulong_direct = 0.5 * e1->innerProduct(e1,q1);
  
#ifdef DEBUG
  double ainv2 = 1./(a*a);
  double erfbeta0over0 = 2*beta/sqrt(MYPI);
  double TWOPI = 2*MYPI;
  int pmax = 10;
  double minusPI2overbeta2 = -MYPI*MYPI/(beta*beta);
  int kmax = 10;
  
  double ulong_direct_constant = 0;
  for (int px = -pmax ; px <= pmax ; px++) {
    for (int py = -pmax ; py <= pmax ; py++) {
      for (int pz = -pmax ; pz <= pmax ; pz++) {
        double r2 = px*px + py*py + pz*pz;
        double r = sqrt(r2);
        double gama = msm4g_smoothing_gama(r/a, nu);
        double erfvalue = r2 == 0 ? erfbeta0over0 : erf(beta*r)/r ;
        ulong_direct_constant +=  (gama/a - erfvalue);
      }
    }
  }
  
  double ulong_fourier_constant = 0;
  for (int kx = -kmax ; kx <= kmax ; kx++) {
    for (int ky = -kmax ; ky <= kmax ; ky++) {
      for (int kz = -kmax ; kz <= kmax ; kz++) {
        if (kx == 0 && ky == 0 && kz == 0) continue;
        double kinvx = kx / Ax ;
        double kinvy = ky / Ay ;
        double kinvz = kz / Az ;
        double k2 = kinvx*kinvx + kinvy*kinvy + kinvz * kinvz;
        double chi = exp(minusPI2overbeta2*k2)/k2 ;
        ulong_fourier_constant += chi ;
      }
    }
  }
  for (int i = 0 ; i < N ; i++) {
    Particle *pi = &(particles[i]);
    pi->acc_long_direct_true[0] = 0;
    pi->acc_long_direct_true[1] = 0;
    pi->acc_long_direct_true[2] = 0;
    pi->acc_long_fourier_true[0] = 0;
    pi->acc_long_fourier_true[1] = 0;
    pi->acc_long_fourier_true[2] = 0;
  }
  so->ulong_direct_true  = 0;
  double complex sum = 0;
  for (int i = 0 ; i < N ; i++) {
    //printf("%d\n",i);
    Particle *pi = &(particles[i]);
    double rix = pi->r.value[0];
    double riy = pi->r.value[1];
    double riz = pi->r.value[2];
    double qi = pi->m ;
    double ulong_direct_sum = 0;
    double complex ulong_fourier_sum = 0;
    for (int j = i+1 ; j < N ; j++) {
      Particle *pj = &(particles[j]);

      double rjx = pj->r.value[0];
      double rjy = pj->r.value[1];
      double rjz = pj->r.value[2];
      double qj = pj->m ;
      for (int px = -pmax ; px <= pmax ; px++) {
        for (int py = -pmax ; py <= pmax ; py++) {
          for (int pz = -pmax ; pz <= pmax ; pz++) {
            double rx = rix - rjx - Ax * px ;
            double ry = riy - rjy - Ay * py ;
            double rz = riz - rjz - Az * pz ;
            double r2 = rx*rx + ry*ry + rz*rz;
            double r = sqrt(r2);
            double gama = msm4g_smoothing_gama(r/a, nu);
            double erfvalue = r2 == 0 ? erfbeta0over0 : erf(beta*r)/r ;
            ulong_direct_sum += qj * (gama/a - erfvalue);
            
            // force calculation
            if (r2 != 0) {
              double magn = qj*(msm4g_smoothing_gamaprime(r/a,nu)*ainv2*r +
                             erfvalue - erfbeta0over0*exp(-beta*beta*r2))/r2;
              pi->acc_long_direct_true[0] += rx * magn ;
              pi->acc_long_direct_true[1] += ry * magn ;
              pi->acc_long_direct_true[2] += rz * magn ;
              pj->acc_long_direct_true[0] -= rx * magn ;
              pj->acc_long_direct_true[1] -= ry * magn ;
              pj->acc_long_direct_true[2] -= rz * magn ;
            }
          }
        }
      }
      for (int kx = -kmax ; kx <= kmax ; kx++) {
        for (int ky = -kmax ; ky <= kmax ; ky++) {
          for (int kz = -kmax ; kz <= kmax ; kz++) {
            if (kx == 0 && ky == 0 && kz == 0) continue;
            double kinvx = kx / Ax ;
            double kinvy = ky / Ay ;
            double kinvz = kz / Az ;
            double dotprod = kinvx * (rix-rjx) + kinvy * (riy-rjy) + kinvz * (riz-rjz);
            double k2 = kinvx*kinvx + kinvy*kinvy + kinvz * kinvz;
            double chi = exp(minusPI2overbeta2*k2)/k2 ;
            double complex value =  qj * chi * cexp(TWOPI * I * dotprod);
            ulong_fourier_sum += value;
            double magn = creal(2.0 * I * value / detA) ;
            pi->acc_long_fourier_true[0] += kinvx * magn ;
            pi->acc_long_fourier_true[1] += kinvy * magn ;
            pi->acc_long_fourier_true[2] += kinvz * magn ;
            pj->acc_long_fourier_true[0] -= kinvx * magn ;
            pj->acc_long_fourier_true[1] -= kinvy * magn ;
            pj->acc_long_fourier_true[2] -= kinvz * magn ;
          }
        }
      }
    }
    so->ulong_direct_true += qi* (ulong_direct_sum + 0.5 * qi*ulong_direct_constant);
    sum += qi * (ulong_fourier_sum + 0.5 * qi*ulong_fourier_constant);
  }

  so->ulong_fourier_true = creal(sum) / (MYPI*detA);
#endif

  AbstractGrid *eLp1 = simulation->gridpotential[L];
  AbstractGrid *qLp1 = simulation->gridmass[L];
  double ulong_fourier = 0.5 * eLp1->innerProduct(eLp1,qLp1);

  so->potentialEnergyShortRange = ushort_real + ushort_csr + ushort_self ;
  so->potentialEnergyLongRange = ulong_direct + ulong_self ;

  so->potentialEnergyTotal = so->potentialEnergyShortRange
      + so->potentialEnergyLongRange;

  so->ushort_real = ushort_real;
  so->ushort_self = ushort_self;
  so->ushort_csr = ushort_csr;
  so->ulong_self = ulong_self;
  so->ulong_direct = ulong_direct - ulong_fourier; /*Four. already included*/
  so->ulong_fourier = ulong_fourier;

  return so->potentialEnergyTotal ;
}

void msm4g_stencil(Simulation *simulation, int l) {
  Boolean periodic = simulation->parameters->periodic ;
  /** @todo Make necessary changes for vacuum boundary. For now
   * it is tailor-made for periodic */
  if (periodic != true) return ;

  int L = simulation->parameters->L ;
  AbstractGrid *stencil = simulation->stencil[l-1];
  int Mx = stencil->nx ;
  int My = stencil->ny ;
  int Mz = stencil->nz ;
  double hx = stencil->hx ;
  double hy = stencil->hy ;
  double hz = stencil->hz ;
  double Ax = simulation->box->wx ;
  double Ay = simulation->box->wy ;
  double Az = simulation->box->wz ;

  int mu = simulation->parameters->mu ;
  int nu = simulation->parameters->nu ;
  double a = simulation->parameters->a ;
  double beta = simulation->parameters->beta ;
  double kmax = simulation->parameters->kmax ;
  double *wprime = simulation->parameters->wprime ;
  //printf("stencil calculation at level:%d\n",l);
  if (l <= L) {
    double al = pow(2,l)*a;

    int Mxmin =   - Mx / 2;
    int Mxmax =    (Mx - 1) / 2;
    int Mymin =   - My / 2;
    int Mymax =    (My - 1) / 2;
    int Mzmin =   - Mz / 2;
    int Mzmax =    (Mz - 1) / 2;

    double *kappa = (double *)calloc(Mx*My*Mz,sizeof(double));
    double *kappahat = (double *)calloc(Mx*My*Mz,sizeof(double));

    for (int nx = Mxmin ; nx <= Mxmax ; nx++) {
      int ni = (nx - Mxmin) * My * Mz ;
      for (int ny = Mymin ; ny <= Mymax ; ny++) {
        int nj = ni + (ny - Mymin) * Mz ;
        for (int nz = Mzmin ; nz <= Mzmax ; nz++) {
          int nk = nj + (nz - Mzmin) ;
          double psum = 0.0;
          int pxmin = - al/Ax + (double)nx / (double) Mx ;
          int pymin = - al/Ay + (double)ny / (double) My ;
          int pzmin = - al/Az + (double)nz / (double) Mz ;
          int pxmax =   al/Ax + (double)nx / (double) Mx ;
          int pymax =   al/Ay + (double)ny / (double) My ;
          int pzmax =   al/Az + (double)nz / (double) Mz ;
          for (int px = pxmin ; px <= pxmax; px++) {
            for (int py = pymin ; py <= pymax; py++) {
              for (int pz = pzmin; pz <= pzmax; pz++) {
                double rx = hx * nx - Ax * px;
                double ry = hy * ny - Ay * py;
                double rz = hz * nz - Az * pz;
                double rlen2 = rx * rx + ry * ry + rz * rz;
                double rlen = sqrt(rlen2);
                double rho = rlen/al ;
                if (rho >= 1) continue;
                double kernel = msm4g_kernel(l,L,rlen,a,beta,nu);
                psum += kernel;
              }
            }
          }
          kappa[nk] = psum ;
          //printf("kappa1 %2d %2d %2d %25.16e\n",nx,ny,nz,psum);
        }
      }
    }
    for (int mx = Mxmin; mx <= Mxmax; mx++) {
      if (abs(mx) > simulation->parameters->nbarx) continue;
      for (int my = Mxmin; my <= Mymax; my++) {
        if (abs(my) > simulation->parameters->nbary) continue;
        for (int mz = Mzmin; mz <= Mzmax; mz++) {
          if (abs(mz) > simulation->parameters->nbarz) continue;
          int kappahatoffset = (mx - Mxmin)* My * Mz  + (my - Mymin) * Mz +
              (mz - Mzmin);
          double sum = 0.0;
          for (int kx = - mu - nu / 2; kx <= mu + nu / 2; kx++) {
            for (int ky = - mu - nu / 2; ky <= mu + nu / 2; ky++) {
              for (int kz = - mu - nu / 2; kz <= mu + nu / 2; kz++) {
                double omega =
                    wprime[abs(kx)] * wprime[abs(ky)] * wprime[abs(kz)];

                int nx = mx + kx ;
                int ny = my + ky ;
                int nz = mz + kz ;
                if (nx < Mxmin) do { nx += Mx; } while (nx < Mxmin);
                if (ny < Mymin) do { ny += My; } while (ny < Mymin);
                if (nz < Mzmin) do { nz += Mz; } while (nz < Mzmin);
                if (nx > Mxmax) do { nx -= Mx; } while (nx > Mxmax);
                if (ny > Mymax) do { ny -= My; } while (ny > Mymax);
                if (nz > Mzmax) do { nz -= Mz; } while (nz > Mzmax);
                int kappaoffset = (nx - Mxmin)* My * Mz  + (ny - Mymin) * Mz +
                    (nz - Mzmin);
                sum += omega * kappa[kappaoffset];
              }
            }
          }
          kappahat[kappahatoffset] = sum ;
          stencil->setElement(stencil,mx-Mxmin,my-Mymin,mz-Mzmin,sum);
        }
      }
    }
    free(kappa);
    free(kappahat);
  } else if (l == L + 1) {
    double detA = Ax * Ay * Az ;
    int Mxmin =   - Mx / 2;
    int Mxmax =    (Mx - 1) / 2;
    int Mymin =   - My / 2;
    int Mymax =    (My - 1) / 2;
    int Mzmin =   - Mz / 2;
    int Mzmax =    (Mz - 1) / 2;
    int cmax = ceil(MSM4G_MAX3(Ax, Ay, Az) * kmax);
    double *cvec = (double *)calloc(2*cmax+1,sizeof(double));
    for (int cn = -cmax ; cn <= cmax ; cn++) {
      cvec[cn+cmax] = msm4g_util_calculate_c(cn, Mx, nu,mu,wprime);
    }
    for (int mx = Mxmin; mx <= Mxmax; mx++) {
      for (int my = Mymin; my <= Mymax; my++) {
        for (int mz = Mzmin; mz <= Mzmax; mz++) {
          double sum = 0.0;
          for (int kx = - Ax * kmax; kx <= Ax*kmax; kx++) {
            for (int ky = - Ay * kmax; ky <= Ay * kmax; ky++) {
              for (int kz = - Az * kmax; kz <= Az * kmax; kz++) {
                if (kx == 0 && ky == 0 && kz == 0)
                  continue;
                double kvecx = kx / Ax; // kvec = inv(A) * k
                double kvecy = ky / Ay;
                double kvecz = kz / Az;
                double k2 = kvecx * kvecx + kvecy * kvecy + kvecz * kvecz;
                // Eq. 5
                double chi = (1.0 / (MYPI * k2 * detA))
                * exp(-MYPI * MYPI * k2 / (beta * beta));
                double dotprod = kx * mx / (double) Mx + ky * my / (double) My
                    + kz * mz / (double) Mz;

                double cx = cvec[kx + cmax];
                double cy = cvec[ky + cmax];
                double cz = cvec[kz + cmax];

                sum += chi * cx * cy * cz * cos(2 * MYPI * dotprod);
              }
            }
          }
          stencil->setElement(stencil,mx-Mxmin,my-Mymin,mz-Mzmin,sum);
        }
      }
    }
    free(cvec);
  }
}

Simulation *msm4g_simulation_new(char *datafile,SimulationBox *box,
    Boolean periodic,int order,double abar,int mu,int L,int Mx,int My,int Mz) {
  Simulation *simulation;
  SimulationParameters *sp;
  SimulationOutput *output;
  Particle *particles;
  int N = 0;

  simulation = (Simulation *)calloc(1,sizeof(Simulation));
  sp = (SimulationParameters *)calloc(1,sizeof(SimulationParameters));
  output = (SimulationOutput *)calloc(1,sizeof(SimulationOutput));
  particles = msm4g_particle_read(datafile,&N);
  double vol = box->wx * box->wy * box->wz ;
  sp->N = N;
  sp->abar = abar;
  sp->mu = mu;
  sp->periodic = periodic;
  sp->nu = order ;
  if ( L== 0 )
    L = floor(1 + 0.5 * log(sp->N)/log(8));
  sp->L = L ;

  /* Number of grid on finest level */
  if (Mx == 0)
    Mx = pow(2,sp->L-1) * ceil(pow(2,1-sp->L)*box->wx * pow(sp->N/(vol),1/3.0));
  if (My == 0)
    My = pow(2,sp->L-1) * ceil(pow(2,1-sp->L)*box->wy * pow(sp->N/(vol),1/3.0));
  if (Mz == 0)
    Mz = pow(2,sp->L-1) * ceil(pow(2,1-sp->L)*box->wz * pow(sp->N/(vol),1/3.0));

  sp->Mx = Mx ;
  sp->My = My ;
  sp->Mz = Mz ;

  /** @todo These are slightly different from the manuscript which uses
   * Mxmin =   -M/2, Mxmax = M/2 - 1 */
  sp->Mxmin =  - Mx / 2;
  sp->Mxmax =   (Mx - 1)/2 ;
  sp->Mymin =  - My / 2;
  sp->Mymax =   (My - 1)/2 ;
  sp->Mzmin =  - Mz / 2;
  sp->Mzmax =   (Mz - 1)/2 ;

  sp->hx = box->wx / sp->Mx ;
  sp->hy = box->wy / sp->My ;
  sp->hz = box->wz / sp->Mz ;
  sp->h  = MSM4G_MAX3(sp->hx,sp->hy,sp->hz);
  sp->a = sp->abar * sp->h ;

  sp->nbarx = ceil(2.0 * sp->a / sp->hx );
  sp->nbary = ceil(2.0 * sp->a / sp->hy );
  sp->nbarz = ceil(2.0 * sp->a / sp->hz );

  double aL = pow(2,L) * sp->a ;
  double detA = box->wx * box->wy * box->wz ;
  double h0 = pow(detA/N,1./3.);
  sp->beta = msm4g_util_choose_beta(aL,TOL_DIRECT,h0) ;
  sp->kmax = msm4g_util_choose_kmax(sp->beta, TOL_FOURIER, h0);
  {
    int l,Mx,My,Mz;
    double hx,hy,hz;
    int extension = sp->nu - 1;
    if (periodic == true)
      extension = 0;

    for (l = 1 ; l <= sp->L ; l++) {
      int Mx = sp->Mx / pow(2,l-1) ;
      int My = sp->My / pow(2,l-1) ;
      int Mz = sp->Mz / pow(2,l-1) ;
      double hx = sp->hx * pow(2,l-1);
      double hy = sp->hy * pow(2,l-1);
      double hz = sp->hz * pow(2,l-1);
//      int nbarx = MSM4G_MIN2(sp->nbarx , Mx / 2);
//      int nbary = MSM4G_MIN2(sp->nbary , My / 2);
//      int nbarz = MSM4G_MIN2(sp->nbarz , Mz / 2);
      simulation->stencil[l-1] = msm4g_grid_dense_new(Mx + extension,
          My + extension,Mz + extension,hx,hy,hz);
      simulation->gridpotential[l-1] = msm4g_grid_dense_new(Mx+extension,
          My+extension,Mz+extension,hx,hy,hz);
      simulation->gridmass[l-1]      = msm4g_grid_dense_new(Mx+extension,
          My+extension,Mz+extension,hx,hy,hz);
    }
    Mx = sp->Mx / pow(2,L-1) ;
    My = sp->My / pow(2,L-1) ;
    Mz = sp->Mz / pow(2,L-1) ;
    hx = sp->hx * pow(2,L-1) ;
    hy = sp->hy * pow(2,L-1) ;
    hz = sp->hz * pow(2,L-1) ;
    simulation->stencil[L] = msm4g_grid_dense_new(Mx+extension,
        My+extension,Mz+extension,hx,hy,hz);
    simulation->gridpotential[L] = msm4g_grid_dense_new(Mx+extension,
        My+extension,Mz+extension,hx,hy,hz);
    simulation->gridmass[L]      = msm4g_grid_dense_new(Mx+extension,
        My+extension,Mz+extension,hx,hy,hz);
  }
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
  LinkedList *binlist = msm4g_bin_generate(box,particles,
      simulation->parameters->N,sp->a,sp->periodic);
  int L = simulation->parameters->L ;
  int N = simulation->parameters->N ;
  msm4g_tic();
  sp->wprime = msm4g_util_omegaprime(sp->mu,sp->nu);
  simulation->output->time_omegaprime = msm4g_toc();
  msm4g_tic();
  msm4g_force_short(binlist, sp->a, simulation);
  simulation->output->time_short_direct = msm4g_toc();
  msm4g_tic();
  msm4g_anterpolation(simulation);
  simulation->output->time_anterpolation = msm4g_toc();

  msm4g_tic();
  for (int l = 2 ; l <= L ; l++) {
    msm4g_restriction(simulation,l);
  }
  simulation->output->time_restriction = msm4g_toc();

  AbstractGrid *qL = simulation->gridmass[L-1];
  AbstractGrid *qLplus1 = simulation->gridmass[L];
  qLplus1->add(qLplus1,qL);

  msm4g_tic();
  for (int l = 1 ; l <= sp->L; l++) {
    msm4g_stencil(simulation,l);
  }
  simulation->output->time_stencil = msm4g_toc();

  msm4g_tic();
  msm4g_stencil(simulation,sp->L + 1);
  simulation->output->time_stencil_fourier = msm4g_toc();


  /* Grid-to-grid mapping */
  for (int l = 1 ; l <= sp->L + 1; l++) {
    AbstractGrid *stencil = simulation->stencil[l-1];
    AbstractGrid *gridmass = simulation->gridmass[l-1];
    AbstractGrid *gridpotential = simulation->gridpotential[l-1];
    msm4g_tic();
    msm4g_grid_potential(stencil,gridmass,gridpotential);
    simulation->output->time_grid_to_grid[l-1] = msm4g_toc();
  }
  /*
  printf("Grid mass\n");
  msm4g_grid_print(simulation->gridmass[0]);
  printf("Stencil\n");
  msm4g_grid_print(simulation->stencil[0]);
  printf("Stencil at 1\n");
    msm4g_grid_print(simulation->stencil[1]);
  printf("Grid potential\n");
  msm4g_grid_print(simulation->gridpotential[0]);
  printf("Energy\n");
  printf("%25.16e\n",0.5*simulation->gridpotential[0]->innerProduct(
      simulation->gridpotential[0],simulation->gridmass[0])); */
  /* Prolongation */
  msm4g_tic();
  msm4g_prolongation(simulation);
  simulation->output->time_prolongation = msm4g_toc();

  /* Interpolation */
  msm4g_tic();
  msm4g_interpolation(simulation);
  simulation->output->time_interpolation = msm4g_toc();

  msm4g_tic();
  msm4g_potential_energy(simulation);
  simulation->output->time_energy = msm4g_toc();


  for (int i = 0 ; i < N ; i++) {
    double shortx = simulation->particles[i].acc_short[0];
    double shorty = simulation->particles[i].acc_short[1];
    double shortz = simulation->particles[i].acc_short[2];
    double longx  = simulation->particles[i].acc_long[0];
    double longy  = simulation->particles[i].acc_long[1];
    double longz  = simulation->particles[i].acc_long[2];
    simulation->particles[i].acc_total[0] = shortx + longx ;
    simulation->particles[i].acc_total[1] = shorty + longy ;
    simulation->particles[i].acc_total[2] = shortz + longz ;
  }

  msm4g_bin_destroy(binlist);
}

void msm4g_simulation_save(Simulation *simulation,FILE *fp) {
  SimulationParameters *sp = simulation->parameters;
  SimulationBox *box = simulation->box;
  SimulationOutput *so = simulation->output;
  int i,l;
  fprintf(fp,"simulation:\n");
  fprintf(fp,"  parameters:\n");
  fprintf(fp,"    abar: %25.16e\n",sp->abar);
  fprintf(fp,"    nu: %d\n",sp->nu);
  fprintf(fp,"    mu: %d\n",sp->mu);
  fprintf(fp,"    periodic: %d\n",sp->periodic);
  fprintf(fp,"    a: %25.16e\n",sp->a);
  fprintf(fp,"    h: %25.16e\n",sp->h);
  fprintf(fp,"    hx: %25.16e\n",sp->hx);
  fprintf(fp,"    hy: %25.16e\n",sp->hy);
  fprintf(fp,"    hz: %25.16e\n",sp->hz);
  fprintf(fp,"    Mx: %d\n",sp->Mx);
  fprintf(fp,"    My: %d\n",sp->My);
  fprintf(fp,"    Mz: %d\n",sp->Mz);
  fprintf(fp,"    nbarx: %d\n",sp->nbarx);
  fprintf(fp,"    nbary: %d\n",sp->nbary);
  fprintf(fp,"    nbarz: %d\n",sp->nbarz);
  fprintf(fp,"    L: %d\n",sp->L);
  fprintf(fp,"    N: %d\n",sp->N);
  fprintf(fp,"    Mxmin: %d\n",sp->Mxmin);
  fprintf(fp,"    Mxmax: %d\n",sp->Mxmax);
  fprintf(fp,"    Mymin: %d\n",sp->Mymin);
  fprintf(fp,"    Mymax: %d\n",sp->Mymax);
  fprintf(fp,"    Mzmin: %d\n",sp->Mzmin);
  fprintf(fp,"    Mzmax: %d\n",sp->Mzmax);
  fprintf(fp,"    wprime: [");
  for (i = 0 ; i < sp->mu + (sp->nu)/2 ; i++)
    fprintf(fp,"%25.16e, ",sp->wprime[i]);
  fprintf(fp,"%25.16e] \n",sp->wprime[i]);
  double sumwprime = sp->wprime[0];
  for (i = 1 ; i < sp->mu + (sp->nu)/2 + 1; i++)
    sumwprime += 2 * sp->wprime[i];
  fprintf(fp,"    sumwprime: %25.16f\n",sumwprime);
  fprintf(fp,"    beta: %25.16lf\n",sp->beta);
  fprintf(fp,"    kmax: %25.16lf\n",sp->kmax);
  fprintf(fp,"  box:\n");
  fprintf(fp,"    x: %25.16e\n",box->x);
  fprintf(fp,"    y: %25.16e\n",box->y);
  fprintf(fp,"    z: %25.16e\n",box->z);
  fprintf(fp,"    wx: %25.16e\n",box->wx);
  fprintf(fp,"    wy: %25.16e\n",box->wy);
  fprintf(fp,"    wz: %25.16e\n",box->wz);
  fprintf(fp,"    stencilsizes: [");
  for (l = 0 ; l < simulation->parameters->L ; l++)
    fprintf(fp,"%d,",simulation->stencil[l]->nx);
  fprintf(fp,"%d]\n",simulation->stencil[l]->nx);


  fprintf(fp,"  output:\n");
  fprintf(fp,"    potentialEnergyShortRange: %25.16e\n",
      so->potentialEnergyShortRange);
  fprintf(fp,"    potentialEnergyLongRange: %25.16e\n",
      so->potentialEnergyLongRange);
  fprintf(fp,"    uconstant: %25.16e\n",so->ushort_csr+so->ushort_self+
      so->ulong_self);
  fprintf(fp,"    utotal: %25.16e\n",so->potentialEnergyTotal);
  fprintf(fp,"    ushort_real: %25.16e\n",so->ushort_real);
#ifdef DEBUG
  fprintf(fp,"    ushort_real_true: %25.16e\n",so->ushort_real_true);
  fprintf(fp,"    ushort_real_error: %25.16e\n",
          fabs(so->ushort_real-so->ushort_real_true)/fabs(so->ushort_real_true));
#endif
  fprintf(fp,"    ushort_self: %25.16e\n",so->ushort_self);
  fprintf(fp,"    ushort_csr: %25.16e\n",so->ushort_csr);
  fprintf(fp,"    ulong_self: %25.16e\n",so->ulong_self);
  fprintf(fp,"    ulong_direct: %25.16e\n",so->ulong_direct);
#ifdef DEBUG
  fprintf(fp,"    ulong_direct_true: %25.16e\n",so->ulong_direct_true);
  fprintf(fp,"    ulong_direct_error: %25.16e\n",
          fabs(so->ulong_direct-so->ulong_direct_true)/fabs(so->ulong_direct_true));
#endif
  fprintf(fp,"    ulong_fourier: %25.16e\n",so->ulong_fourier);
#ifdef DEBUG
  int N = simulation->parameters->N;

  fprintf(fp,"    ulong_fourier_true: %25.16e\n",so->ulong_fourier_true);
  fprintf(fp,"    ulong_fourier_error: %25.16e\n",
          fabs(so->ulong_fourier-so->ulong_fourier_true)/fabs(so->ulong_fourier_true));
  
  double avg_force_error = 0;
  for (int i = 0 ; i < N ; i++) {
    Particle *pi = &(simulation->particles[i]);
    double fx_true = pi->acc_short_true[0] + pi->acc_long_direct_true[0] + pi->acc_long_fourier_true[0];
    double fy_true = pi->acc_short_true[1] + pi->acc_long_direct_true[1] + pi->acc_long_fourier_true[1];
    double fz_true = pi->acc_short_true[2] + pi->acc_long_direct_true[2] + pi->acc_long_fourier_true[2];
    double fx = pi->acc_total[0];
    double fy = pi->acc_total[1];
    double fz = pi->acc_total[2];
    double force_error = sqrt((fx_true-fx)*(fx_true-fx)+(fy_true-fy)*(fy_true-fy)+(fz_true-fz)*(fz_true-fz))/
                              sqrt(fx*fx+fy*fy+fz*fz);
    avg_force_error += force_error ;
  }
  avg_force_error /= N ;
  fprintf(fp,"    avg_force_error: %25.16e\n",avg_force_error);
#endif
  fprintf(fp,"    time_omegaprime: %25.16lf\n",so->time_omegaprime);
  fprintf(fp,"    time_anterpolation: %25.16lf\n",so->time_anterpolation);
  fprintf(fp,"    time_interpolation: %25.16lf\n",so->time_interpolation);
  fprintf(fp,"    time_short_direct: %25.16lf\n",so->time_short_direct);
  fprintf(fp,"    time_restriction: %25.16lf\n",so->time_restriction);
  fprintf(fp,"    time_prolongation: %25.16lf\n",so->time_prolongation);
  fprintf(fp,"    time_energy: %25.16lf\n",so->time_energy);
  fprintf(fp,"    time_stencil: %25.16lf\n",so->time_stencil);
  fprintf(fp,"    time_stencil_fourier: %25.16lf\n",so->time_stencil_fourier);
  fprintf(fp,"    time_grid_to_grid: [");
  for (l = 1 ; l < simulation->parameters->L ; l++)
    fprintf(fp,"%25.16lf,",so->time_grid_to_grid[l-1]);
  fprintf(fp,"%25.16lf]\n",so->time_grid_to_grid[l-1]);
  fprintf(fp,"    shortRangeInteractionCount: %d\n",
      so->shortRangeInteractionCount);
  fprintf(fp,"    kernelHatEvaluationsNeeded: [");
  for (l = 1 ; l < simulation->parameters->L ; l++)
    fprintf(fp,"%d,",so->kernelHatEvaluationsNeeded[l-1]);
  fprintf(fp,"%d]\n",so->kernelHatEvaluationsNeeded[l-1]);
  fprintf(fp,"    kernelHatEvaluationsComputed: [");
  for (l = 1 ; l < simulation->parameters->L ; l++)
    fprintf(fp,"%d,",so->kernelHatEvaluationsComputed[l-1]);
  fprintf(fp,"%d]\n",so->kernelHatEvaluationsComputed[l-1]);

  fprintf(fp,"    kernelEvaluationsNeeded: [");
  for (l = 1 ; l < simulation->parameters->L ; l++)
    fprintf(fp,"%d,",so->kernelEvaluationsNeeded[l-1]);
  fprintf(fp,"%d]\n",so->kernelEvaluationsNeeded[l-1]);
  fprintf(fp,"    kernelEvaluationsComputed: [");
  for (l = 1 ; l < simulation->parameters->L ; l++)
    fprintf(fp,"%d,",so->kernelEvaluationsComputed[l-1]);
  fprintf(fp,"%d]\n",so->kernelEvaluationsComputed[l-1]);
}

void msm4g_simulation_delete(Simulation *simulation) {
  for (int l = 0 ; l <= simulation->parameters->L ; l++) {
    if (simulation->gridmass[l] != NULL)
      msm4g_grid_destroy(&(simulation->gridmass[l]));
    if (simulation->gridpotential[l] != NULL)
      msm4g_grid_destroy(&(simulation->gridpotential[l]));
    if (simulation->stencil[l] != NULL)
      msm4g_grid_destroy(&(simulation->stencil[l]));
  }

  msm4g_box_destroy(simulation->box);
  free(simulation->parameters->wprime);
  free(simulation->particles);
  free(simulation->parameters);
  free(simulation->output);
  free(simulation);
}

void msm4g_anterpolation(Simulation *simulation) {
  int i0,j0,k0;
  int p = simulation->parameters->nu ;
  double x,y,z;
  double x0,y0,z0;
  double rx_hx,ry_hy,rz_hz;
  int s_edge;
  double h = simulation->parameters->h ;
  int Mx = simulation->parameters->Mx ;
  int My = simulation->parameters->My ;
  int Mz = simulation->parameters->Mz ;

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

  grid = simulation->gridmass[0];
  if (simulation->parameters->periodic == true) {
    double hx = simulation->parameters->hx ;
    double hy = simulation->parameters->hy ;
    double hz = simulation->parameters->hz ;
    for (int i = 0 ; i < N ; i++) {
      double rx = particles[i].r.value[0];
      double ry = particles[i].r.value[1];
      double rz = particles[i].r.value[2];
      double mass  = particles[i].m ;
      double mx = floor(rx/hx) ;
      double my = floor(ry/hy) ;
      double mz = floor(rz/hz) ;
      double tx = rx/hx - mx ;
      double ty = ry/hy - my ;
      double tz = rz/hz - mz ;
      for (int nx = 1 - nu/2 ; nx <= nu/2 ; nx++) {
        double phix = msm4g_bases_bspline(nu, tx - nx + nu / 2);
        int mnx = mx + nx ;
        if (mnx <   0) do { mnx += Mx; } while (mnx <   0);
        if (mnx >= Mx) do { mnx -= Mx; } while (mnx >= Mx);
        for (int ny = 1 - nu/2 ; ny <= nu/2 ; ny++) {
          double phiy = msm4g_bases_bspline(nu, ty - ny + nu / 2);
          int mny = my + ny ;
          if (mny <   0) do { mny += My; } while (mny <   0);
          if (mny >= My) do { mny -= My; } while (mny >= My);
          for (int nz = 1 - nu/2 ; nz <= nu/2 ; nz++) {
            double phiz = msm4g_bases_bspline(nu, tz - nz + nu / 2);
            int mnz = mz + nz ;
            if (mnz <   0) do { mnz += Mz; } while (mnz <   0);
            if (mnz >= Mz) do { mnz -= Mz; } while (mnz >= Mz);
            double oldValue = grid->getElement(grid,mnx,mny,mnz);
            double newValue = oldValue + mass * phix * phiy * phiz ;
            grid->setElement(grid,mnx,mny,mnz,newValue);
          }
        }
      }
    }
  } else {
    s_edge = p/2-1;
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

      for (int v = 0; v < nu ; v++) {
        tx = rx_hx - i0 - v;
        ty = ry_hy - j0 - v;
        tz = rz_hz - k0 - v;
        phix[v] = msm4g_bases_bspline(nu,tx+nu/2) ;
        phiy[v] = msm4g_bases_bspline(nu,ty+nu/2) ;
        phiz[v] = msm4g_bases_bspline(nu,tz+nu/2) ;
      }

      for (int i = 0; i < nu ; i++) {
        for (int j = 0 ; j < nu ; j++) {
          for (int k = 0; k < nu ; k++) {
            double g = phix[i] * phiy[j] * phiz[k] * mass;
            double gold = grid->getElement(grid,i+i0,j+j0,k+k0);
            grid->setElement(grid,i+i0,j+j0,k+k0,g+gold);
          }
        }
      }
    }
  }
}

void msm4g_restriction(Simulation *simulation, int l) {
  AbstractGrid *qm = simulation->gridmass[l-1];  /* Coarse level */
  AbstractGrid *qn = simulation->gridmass[l-2];/* Fine level */
  int nu = simulation->parameters->nu ;
  int Mx = qm->nx ;
  int My = qm->ny ;
  int Mz = qm->nz ;
  int Nx = qn->nx ;
  int Ny = qn->ny ;
  int Nz = qn->nz ;
  for (int mx = 0 ; mx < Mx ; mx++) {
    for (int my = 0 ; my < My ; my++) {
      for (int mz = 0 ; mz < Mz ; mz++) {
        double sum = 0.0;
        int nxmin = 2 * mx - nu / 2 ;
        int nymin = 2 * my - nu / 2 ;
        int nzmin = 2 * mz - nu / 2 ;
        int nxmax = 2 * mx + nu / 2 ;
        int nymax = 2 * my + nu / 2 ;
        int nzmax = 2 * mz + nu / 2 ;
        for (int nx = nxmin ; nx <= nxmax ; nx++) {
          int nxwrapped = nx ;
          if (nxwrapped <  0 ) do {nxwrapped += Nx;} while (nxwrapped <  0 );
          if (nxwrapped >= Nx) do {nxwrapped -= Nx;} while (nxwrapped >= Nx);
          double jnx = msm4g_util_jn(nu, nx-2*mx) ;
          for (int ny = nymin ; ny <= nymax ; ny++) {
            int nywrapped = ny ;
            if (nywrapped <  0 ) do {nywrapped += Ny;} while (nywrapped <  0 );
            if (nywrapped >= Ny) do {nywrapped -= Ny;} while (nywrapped >= Ny);
            double jny = msm4g_util_jn(nu, ny-2*my) ;
            for (int nz = nzmin ; nz <= nzmax ; nz++) {
              int nzwrapped = nz ;
              if (nzwrapped<  0 ) do {nzwrapped += Nz;} while (nzwrapped <  0 );
              if (nzwrapped>= Nz) do {nzwrapped -= Nz;} while (nzwrapped >= Nz);
              double jnz = msm4g_util_jn(nu, nz-2*mz) ;
              double jn = jnx * jny * jnz ;
              sum += jn * qn->getElement(qn,nxwrapped,nywrapped,nzwrapped);
            }
          }
        }
        qm->setElement(qm,mx,my,mz,sum);
      }
    }
  }
}

void msm4g_prolongation(Simulation *simulation) {
  int L = simulation->parameters->L ;
  int nu = simulation->parameters->nu ;
  AbstractGrid *eLplus1 = simulation->gridpotential[L];
  AbstractGrid *eL      = simulation->gridpotential[L-1];
  /* eL = eL + eLplus1 */
  eL->add(eL,eLplus1);

  for (int l = L-1 ; l >= 1 ; l--) {
    AbstractGrid *en = simulation->gridpotential[l];/* Coarse level */
    AbstractGrid *em = simulation->gridpotential[l-1];  /* Fine level */
    int Mx = em->nx ;
    int My = em->ny ;
    int Mz = em->nz ;
    int Nx = en->nx ;
    int Ny = en->ny ;
    int Nz = en->nz ;
    for (int mx = 0 ; mx < Mx ; mx++) {
      for (int my = 0 ; my < My ; my++) {
        for (int mz = 0 ; mz < Mz ; mz++) {
          double sum = 0.0;
          int nxmin = mx / 2.0 - nu / 4.0 ;
          int nymin = my / 2.0 - nu / 4.0 ;
          int nzmin = mz / 2.0 - nu / 4.0 ;
          int nxmax = mx / 2.0 + nu / 4.0 ;
          int nymax = my / 2.0 + nu / 4.0 ;
          int nzmax = mz / 2.0 + nu / 4.0 ;
          for (int nx = nxmin ; nx <= nxmax ; nx++) {
            int nxwrapped = nx ;
            if (nxwrapped <  0 ) do {nxwrapped+= Nx;} while (nxwrapped <  0 );
            if (nxwrapped >= Nx) do {nxwrapped-= Nx;} while (nxwrapped >= Nx);
            double jnx = msm4g_util_jn(nu, mx-2*nx) ;
            for (int ny = nymin ; ny <= nymax ; ny++) {
              int nywrapped = ny ;
              if (nywrapped <  0 ) do {nywrapped += Ny;} while (nywrapped<  0 );
              if (nywrapped >= Ny) do {nywrapped -= Ny;} while (nywrapped>= Ny);
              double jny = msm4g_util_jn(nu, my-2*ny) ;
              for (int nz = nzmin ; nz <= nzmax ; nz++) {
                int nzwrapped = nz ;
                if (nzwrapped<  0 ) do {nzwrapped+= Nz;} while (nzwrapped<  0 );
                if (nzwrapped>= Nz) do {nzwrapped-= Nz;} while (nzwrapped>= Nz);
                double jnz = msm4g_util_jn(nu, mz-2*nz) ;
                double jn = jnx * jny * jnz ;
                sum += jn * en->getElement(en,nxwrapped,nywrapped,nzwrapped);
              }
            }
          }
          double oldvalue = em->getElement(em,mx,my,mz);
          double newvalue = oldvalue + sum ;
          em->setElement(em,mx,my,mz,newvalue); /** @todo Clean it */
        }
      }
    }
  }
}

void msm4g_interpolation(Simulation *simulation) {
  AbstractGrid * grid = simulation->gridpotential[0];
  if (simulation->parameters->periodic == true) {
    double hx = simulation->parameters->hx ;
    double hy = simulation->parameters->hy ;
    double hz = simulation->parameters->hz ;
    Particle *particles = simulation->particles;
    int N = simulation->parameters->N ;
    int nu = simulation->parameters->nu ;
    int Mx = simulation->parameters->Mx;
    int My = simulation->parameters->My;
    int Mz = simulation->parameters->Mz;
    for (int i = 0 ; i < N ; i++) {
      double rx = particles[i].r.value[0];
      double ry = particles[i].r.value[1];
      double rz = particles[i].r.value[2];
      double mx = floor(rx/hx) ;
      double my = floor(ry/hy) ;
      double mz = floor(rz/hz) ;
      double tx = rx/hx - mx ;
      double ty = ry/hy - my ;
      double tz = rz/hz - mz ;
      double varphidx = 0.0;
      double varphidy = 0.0;
      double varphidz = 0.0;
      for (int nx = 1 - nu/2 ; nx <= nu/2 ; nx++) {
        double  phix = msm4g_bases_bspline(nu, tx - nx + nu / 2);
        double dphix = msm4g_bases_bsplineprime(nu, tx - nx + nu / 2);
        int mnx = mx + nx ;
        if (mnx <   0) do { mnx += Mx; } while (mnx <   0);
        if (mnx >= Mx) do { mnx -= Mx; } while (mnx >= Mx);
        for (int ny = 1 - nu/2 ; ny <= nu/2 ; ny++) {
          double  phiy = msm4g_bases_bspline(nu, ty - ny + nu / 2);
          double dphiy = msm4g_bases_bsplineprime(nu, ty - ny + nu / 2);
          int mny = my + ny ;
          if (mny <   0) do { mny += My; } while (mny <   0);
          if (mny >= My) do { mny -= My; } while (mny >= My);
          for (int nz = 1 - nu/2 ; nz <= nu/2 ; nz++) {
            double  phiz = msm4g_bases_bspline(nu, tz - nz + nu / 2);
            double dphiz = msm4g_bases_bsplineprime(nu, tz - nz + nu / 2);
            int mnz = mz + nz ;
            if (mnz <   0) do { mnz += Mz; } while (mnz <   0);
            if (mnz >= Mz) do { mnz -= Mz; } while (mnz >= Mz);
            double gridPotential = grid->getElement(grid,mnx,mny,mnz);
            varphidx +=  gridPotential * dphix *  phiy *  phiz / hx;
            varphidy +=  gridPotential *  phix * dphiy *  phiz / hy;
            varphidz +=  gridPotential *  phix *  phiy * dphiz / hz;
          }
        }
      }
      particles[i].acc_long[0] = varphidx ;
      particles[i].acc_long[1] = varphidy ;
      particles[i].acc_long[2] = varphidz ;
    }
  }
}

void msm4g_grid_print(AbstractGrid *grid) {
  for (int i = 0 ; i < grid->nx ; i++) {
    for (int j = 0 ; j < grid->ny ; j++) {
      for (int k = 0 ; k < grid->nz ; k++) {
        printf("%2d %2d %2d %25.16e\n",i,j,k,grid->getElement(grid,i,j,k));
      }
    }
  }
}

void msm4g_grid_destroy(AbstractGrid **grid) {
  if (grid != NULL) {
    (*grid)->destructor(grid);
    free(*grid);
  }
  *grid = NULL;
}

AbstractGrid *msm4g_grid_dense_new(int nx, int ny, int nz,double hx,
    double hy,double hz) {
  AbstractGrid *grid;
  DenseGrid *densegrid;

  densegrid = malloc(sizeof(DenseGrid));
  grid = (AbstractGrid *)densegrid;
  grid->hx = hx;
  grid->hy = hy;
  grid->hz = hz;
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->constructor = msm4g_grid_dense_new;
  grid->destructor  = msm4g_grid_dense_destroy;
  grid->setElement  = msm4g_grid_dense_setElement;
  grid->getElement  = msm4g_grid_dense_getElement;
  grid->reset       = msm4g_grid_dense_reset;
  grid->innerProduct= msm4g_grid_dense_innerProduct;
  grid->add         = msm4g_grid_dense_add;
  grid->sum         = msm4g_grid_dense_sum;

  densegrid->data = calloc(nx*ny*nz,sizeof(double));

  return grid;
}

void msm4g_grid_dense_setElement(AbstractGrid *grid,int i,int j,
    int k,double value) {
  DenseGrid *densegrid = (DenseGrid *)grid;
  int position;

  position = i * grid->ny * grid->nz + j * grid->nz + k;
  densegrid->data[position] = value;
}

double msm4g_grid_dense_getElement(AbstractGrid *grid,int i,int j,int k) {
  DenseGrid *densegrid = (DenseGrid *) grid;
  int position;

  position = i * grid->ny * grid->nz + j * grid->nz + k;
  return densegrid->data[position];
}

double msm4g_grid_dense_innerProduct(AbstractGrid *grid1,AbstractGrid *grid2){
  int size1 = grid1->nx * grid1->ny * grid1->nz ;
  int size2 = grid2->nx * grid2->ny * grid2->nz ;
  if (size1 != size2) {
    fprintf(stderr,"Size of the grids are not same\n");
    return 0;
  }

  DenseGrid *denseGrid1 = (DenseGrid *)grid1;
  DenseGrid *denseGrid2 = (DenseGrid *)grid2;
  double innerProduct = 0.0;
  for (int i = 0 ; i < size1 ; i++)
    innerProduct += denseGrid1->data[i] * denseGrid2->data[i];
  return innerProduct;
}

double msm4g_grid_dense_sum(AbstractGrid *grid) {
  DenseGrid *denseGrid = (DenseGrid *)grid;
  double sum = 0.0;
  int n = grid->nx * grid->ny * grid->nz ;
  for (int i = 0 ; i < n; i++) {
    sum += denseGrid->data[i] ;
  }
  return sum;
}

void msm4g_grid_dense_add(AbstractGrid *grid1,AbstractGrid *grid2) {
  int size1 = grid1->nx * grid1->ny * grid1->nz ;
  int size2 = grid2->nx * grid2->ny * grid2->nz ;
  if (size1 != size2) {
    fprintf(stderr,"Size of the grids are not same\n");
    return ;
  }

  DenseGrid *denseGrid1 = (DenseGrid *)grid1;
  DenseGrid *denseGrid2 = (DenseGrid *)grid2;
  for (int i = 0 ; i < size1 ; i++)
    denseGrid1->data[i]  +=  denseGrid2->data[i];
}

void msm4g_grid_dense_reset(AbstractGrid *grid,double value) {
  int i,size;
  DenseGrid *self;
  self = (DenseGrid *)grid;
  size = grid->nx * grid->ny * grid->nz ;
  for (i=0; i<size; i++) {
    self->data[i] = value;
  }
}

void msm4g_grid_dense_destroy(AbstractGrid **grid) {
  DenseGrid *densegrid;
  densegrid = (DenseGrid *)(*grid);
  free(densegrid->data);
  densegrid->data = NULL;
}

void msm4g_grid_potential(AbstractGrid *stencil, AbstractGrid *gridmass,
    AbstractGrid *gridpotential) {
  double *stencildata = ((DenseGrid *) stencil)->data;
  double *gridmassdata = ((DenseGrid *) gridmass)->data;
  double *gridpotentialdata = ((DenseGrid *) gridpotential)->data;

  int mx = gridmass->nx ;
  int my = gridmass->ny ;
  int mz = gridmass->nz ;
  int sx = stencil->nx ;
  int sy = stencil->ny ;
  int sz = stencil->nz ;
  int sxmin =   - sx / 2 ;
  int sxmax =    (sx - 1) / 2 ;
  int symin =   - sy / 2 ;
  int symax =    (sy - 1) / 2 ;
  int szmin =   - sz / 2 ;
  int szmax =    (sz - 1) / 2 ;
//  int sxminlogical =   - mx / 2 ;
//  int sxmaxlogical =    (mx - 1) / 2 ;
//  int syminlogical =   - my / 2 ;
//  int symaxlogical =    (my - 1) / 2 ;
//  int szminlogical =   - mz / 2 ;
//  int szmaxlogical =    (mz - 1) / 2 ;
  for (int mi = 0 ; mi < mx ; mi++) {
    for (int mj = 0 ; mj < my ; mj++) {
      for (int mk = 0 ; mk < mz ; mk++) {
        double potentialsum = 0.0;
        for (int ni = 0 ; ni < mx ; ni++) {
          int mni = mi - ni ;
          if (mni < sxmin) do {mni += sx;} while (mni < sxmin);
          if (mni > sxmax) do {mni -= sx;} while (mni > sxmax);
          int si = mni - sxmin ;
          for (int nj = 0 ; nj < my ; nj++) {
            int mnj = mj - nj ;
            if (mnj < symin) do {mnj += sy;} while (mnj < symin);
            if (mnj > symax) do {mnj -= sy;} while (mnj > symax);
            int sj = mnj - symin ;
            for (int nk = 0 ; nk < mz ; nk++) {
              int mnk = mk - nk ;
              if (mnk < szmin) do {mnk += sz;} while (mnk < szmin);
              if (mnk > szmax) do {mnk -= sz;} while (mnk > szmax);
              int sk = mnk - szmin ;

              int sindex = si * sy * sz + sj * sz + sk  ;
              int gridmassindex = ni * my * mz + nj * mz + nk  ;
              double stencilvalue = stencildata[sindex];
              double gridmassvalue = gridmassdata[gridmassindex];
              potentialsum += stencilvalue * gridmassvalue ;
            }
          }
        }
        int gridpotentialindex = mi * my * mz + mj * mz + mk ;
        gridpotentialdata[gridpotentialindex] = potentialsum ;
      }
    }
  }
}
double msm4g_smoothing_gama(double rho,int nu) {
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

double msm4g_smoothing_gamaprime(double rho,int nu) {
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

double msm4g_kernel(int l,int L, double r,double a,double beta,int nu) {
  double out = 0.0;
  if (l == 0 ) {
    out = 1/r - (1/a)*msm4g_smoothing_gama(r/a,nu);
  } else if (l == L + 1 ){
    if (fabs(r) < DBL_EPSILON)
      out = 2.0 * beta / sqrt(MYPI);
    else
      out = erf(beta * r) / r ;
  } else if (l == L) {
    double erfvalue ;
    if (fabs(r) < DBL_EPSILON)
      erfvalue = 2.0 * beta / sqrt(MYPI);
    else
      erfvalue = erf(beta * r) / r ;
    out = (1./(pow(2,l-1)*a)) *
        msm4g_smoothing_gama(r/(pow(2,l-1)*a),nu) - erfvalue ;
  } else {
    out = (1./(pow(2,l-1)*a)) * msm4g_smoothing_gama(r/(pow(2,l-1)*a),nu) -
        (1./(pow(2,l  )*a)) * msm4g_smoothing_gama(r/(pow(2,l  )*a),nu);
  }
  return out;
}

double msm4g_force_short(LinkedList *binlist,double threshold,
    Simulation *simulation) {
  simulation->output->shortRangeInteractionCount = 0;
  double potential = 0;
#ifdef DEBUG
  SimulationParameters *sp = simulation->parameters;
  int N = sp->N ;
  int order = simulation->parameters->nu;
  double a = simulation->parameters->a ;
  double ainv = 1.0/a;
  double ainv2 = ainv * ainv ;
  double Ax = simulation->box->wx;
  double Ay = simulation->box->wy;
  double Az = simulation->box->wz;
  Particle *particles = simulation->particles;
  for (int i = 0 ; i < N ; i++) {
    Particle *particleI = &(particles[i]) ;
    particleI->acc_short_true[0] = 0;
    particleI->acc_short_true[1] = 0;
    particleI->acc_short_true[2] = 0;
  }
  
  for (int i = 0 ; i < N ; i++) {
    Particle *particleI = &(particles[i]) ;
    double rix = particleI->r.value[0] ;
    double riy = particleI->r.value[1] ;
    double riz = particleI->r.value[2] ;
    double mI = particleI->m ;
    for (int j = i + 1 ; j < N ; j++ ) {
      Particle *particleJ = &(particles[j]);
      double rjx = particleJ->r.value[0] ;
      double rjy = particleJ->r.value[1] ;
      double rjz = particleJ->r.value[2] ;
      double mJ = particleJ->m ;
      int pxmin = (rix - rjx - a)/Ax ;
      int pxmax = (rix - rjx + a)/Ax ;
      int pymin = (riy - rjy - a)/Ay ;
      int pymax = (riy - rjy + a)/Ay ;
      int pzmin = (riz - rjz - a)/Az ;
      int pzmax = (riz - rjz + a)/Az ;
      for (int px = pxmin ; px <= pxmax ; px++) {
        double rx = rix - rjx - Ax * px ;
        double r2i =  rx*rx ;
        double accJx = mJ * rx ;
        double accIx = mI * rx ;
        for (int py = pymin ; py <= pymax ; py++) {
          double ry = riy - rjy - Ay * py ;
          double r2j = r2i + ry*ry;
          double accJy = mJ * ry ;
          double accIy = mI * ry ;
          for (int pz = pzmin ; pz <= pzmax ; pz++) {
            double rz = riz - rjz - Az * pz ;
            double accJz = mJ * rz ;
            double accIz = mI * rz ;
            double r2 = r2j + rz*rz;
            double r = sqrt(r2);
            double rinv = 1.0/r;
            double rovera = r * ainv ;
            double g0 =  rinv -  msm4g_smoothing_gama(rovera,order)*ainv;
            particleI->potential_short_real_true += 0.5 * particleJ->m  * g0;
            particleJ->potential_short_real_true += 0.5 * particleI->m  * g0;
            double gamaprime = (-rinv * rinv - ainv2 *
                msm4g_smoothing_gamaprime(rovera, order))*rinv ;
            particleI->acc_short_true[0] +=  accJx * gamaprime ;
            particleI->acc_short_true[1] +=  accJy * gamaprime ;
            particleI->acc_short_true[2] +=  accJz * gamaprime ;
            particleJ->acc_short_true[0] -=  accIx * gamaprime ;
            particleJ->acc_short_true[1] -=  accIy * gamaprime ;
            particleJ->acc_short_true[2] -=  accIz * gamaprime ;
          }
        }
      }
    }
  }
#endif
  LinkedListElement *currBin=binlist->head;
  while (currBin != NULL) {
    Bin *bin = (Bin *)currBin->data;
    if (bin->isGhost == true) {
      currBin = currBin->next;
      continue;
    }

    //printf("within    [%d %d %d]\n",bin->nx,bin->ny,bin->nz);
    potential += msm4g_force_short_withinBin(bin,threshold,simulation);
    //printf("Bin\n");
    //msm4g_bin_print(bin);

    LinkedListElement *element = bin->neighbors->head;
    while (element != NULL) {
      Bin *neighbor = (Bin *)element->data;
      if ((!neighbor->isGhost  && bin->cantorindex <= neighbor->cantorindex) ||
          ( neighbor->isGhost  &&
              bin->cantorindex <= neighbor->ghostOf->cantorindex)) {
        element = element->next ;
        continue;
      }
      if (neighbor->isGhost) {
        /* Check is ghost neigbhor is a ghost of a neighbor */
        Boolean skipNeighbor = false;
        int searchindex = neighbor->ghostOf->cantorindex;
        LinkedListElement *tmp=bin->neighbors->head;
        while (tmp != NULL) {
          Bin *tmpBin = (Bin *)tmp->data;
          if (tmpBin->cantorindex == searchindex) {
            skipNeighbor = true;
            break;
          }
          tmp = tmp->next;
        }
        if (!skipNeighbor)
          potential += msm4g_force_short_betweenBin(bin,neighbor->ghostOf,
              threshold,simulation);
      } else {
        potential += msm4g_force_short_betweenBin(bin,neighbor,threshold,
            simulation);
      }
      element = element->next ;
    }
    currBin = currBin->next;
  }
  return potential;
}

double msm4g_force_short_withinBin(Bin *bin, double threshold,
    Simulation *simulation) {
  double potentialTotal = 0.0;
  double potential;

  LinkedListElement *currI = bin->particles->head;
  while (currI->next != NULL) {
    Particle *particleI = (Particle *)currI->data;
    LinkedListElement *currJ = currI->next;
    while (currJ != NULL) {
      Particle *particleJ = (Particle *)currJ->data;
      potential = msm4g_force_short_particlePair(particleI,particleJ,
          threshold,simulation);
      potentialTotal += potential;
      currJ = currJ->next;
    }
    currI = currI->next;
  }
  return potentialTotal;
}

double msm4g_force_short_betweenBin(Bin *bin, Bin *neighbor,double threshold,
    Simulation *simulation) {
  LinkedList *particlesI = bin->particles;
  LinkedList *particlesJ = neighbor->particles;
  LinkedListElement *currI, *currJ;
  Particle *particleI,*particleJ;
  double potential = 0.0;
  currI = particlesI->head;
  while (currI != NULL) {
    particleI = (Particle *)currI->data;
    currJ = particlesJ->head;
    while (currJ != NULL) {
      particleJ = (Particle *)currJ->data;
      potential += msm4g_force_short_particlePair(particleI,particleJ,
          threshold,simulation);

      currJ = currJ->next;
    }
    currI = currI->next;
  }
  return potential;
}

double msm4g_force_short_particlePair(Particle *particleI, Particle *particleJ,
    double a, Simulation *simulation) {
  simulation->output->shortRangeInteractionCount++;
  D3Vector rij;
  double r,r2;
  double a2 = a * a;
  double ainv = 1.0/a;
  double ainv2 = ainv * ainv ;
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
      double rx = rix - rjx - Ax * px ;
      double r2i =  rx*rx ;
      double accJx = particleJ->m * rx ;
      double accIx = particleI->m * rx ;
      for (int py = pymin ; py <= pymax ; py++) {
        double ry = riy - rjy - Ay * py ;
        double r2j = r2i + ry*ry;
        double accJy = particleJ->m * ry ;
        double accIy = particleI->m * ry ;
        for (int pz = pzmin ; pz <= pzmax ; pz++) {
          double rz = riz - rjz - Az * pz ;
          double accJz = particleJ->m * rz ;
          double accIz = particleI->m * rz ;
          double r2 = r2j + rz*rz;
          double r = sqrt(r2);
          double rinv = 1.0/r;
          double rovera = r * ainv ;
          double g0 =  rinv -  msm4g_smoothing_gama(rovera,order)*ainv;
          particleI->potential_short_real += 0.5 * particleJ->m  * g0;
          particleJ->potential_short_real += 0.5 * particleI->m  * g0;
          double gamaprime = (-rinv * rinv - ainv2 *
              msm4g_smoothing_gamaprime(rovera, order))*rinv ;
          particleI->acc_short[0] +=  accJx * gamaprime ;
          particleI->acc_short[1] +=  accJy * gamaprime ;
          particleI->acc_short[2] +=  accJz * gamaprime ;
          particleJ->acc_short[0] -=  accIx * gamaprime ;
          particleJ->acc_short[1] -=  accIy * gamaprime ;
          particleJ->acc_short[2] -=  accIz * gamaprime ;
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
    if (r2 < a2) {
      smoothedkernel = 1.0/r -  msm4g_smoothing_gama(r/a,order)/a;
      smoothedkernelderivative = 1.0/r2 -
          msm4g_smoothing_gamaprime(r/a,order)/a2;
      coeff = particleI->m * particleJ->m * smoothedkernelderivative/r;
      /* f_j = f_j + coeff*rij */
      msm4g_d3vector_daxpy(&(particleJ->fshort) ,
          &(particleJ->fshort), coeff, &rij);
      /* f_i = f_i - coeff*rij */
      msm4g_d3vector_daxpy(&(particleI->fshort) ,
          &(particleI->fshort), coeff, &rij);
      potential = particleI->m * particleJ->m * smoothedkernel;
      particleI->potential_short_real += 0.5 * particleJ->m  * smoothedkernel;
      particleJ->potential_short_real += 0.5 * particleI->m  * smoothedkernel;
    }
  }

  return potential;
}

void msm4g_d3vector_set(D3Vector *d3vector,double x,double y,double z) {
  d3vector->value[0] = x;
  d3vector->value[1] = y;
  d3vector->value[2] = z;
}

void msm4g_d3vector_daxpy(D3Vector *z,D3Vector *x,double a,D3Vector *y) {
  int i;
  for (i=0; i<3; i++) {
    z->value[i] = x->value[i] + a*y->value[i];
  }
}

double msm4g_d3vector_norm(D3Vector *x) {
  return sqrt(msm4g_d3vector_normsquare(x));
}

double msm4g_d3vector_normsquare(D3Vector *x) {
  double normsquare = 0.0;
  int i;
  for (i=0; i<3; i++) {
    normsquare += x->value[i] * x->value[i];
  }
  return normsquare;
}

void msm4g_d3vector_print(D3Vector *x) {
  int i;
  for (i=0; i<3; i++) {
    printf("%10.3E ",x->value[i]);
  }
}

void msm4g_i3vector_set(I3Vector *i3vector,int x,int y,int z) {
  i3vector->value[0] = x;
  i3vector->value[1] = y;
  i3vector->value[2] = z;
}

void msm4g_i3vector_copy(I3Vector *to, I3Vector from) {
  to->value[0] = from.value[0];
  to->value[1] = from.value[1];
  to->value[2] = from.value[2];
}

Boolean msm4g_i3vector_isequal(I3Vector *x,I3Vector *y) {
  int i;
  if (x == NULL || y == NULL)
    return false;
  for (i=0;i<3;i++) {
    if (x->value[i] != y->value[i])
      return false;
  }
  return true;
}

void msm4g_i3vector_print(I3Vector *x) {
  printf("[%d,%d,%d] ",x->value[0],x->value[1],x->value[2]);
}

int msm4g_linkedlist_size(LinkedList *list) {
  int size = 0;
  LinkedListElement *curr;
  curr = list->head;
  while (curr != NULL) {
    size = size + 1;
    curr=curr->next;
  }
  return size;
}

LinkedList *msm4g_linkedlist_new() {
  LinkedList *newlist ;
  newlist = malloc(sizeof(LinkedList));
  newlist->head = NULL;
  newlist->tail = NULL;
  return newlist;
}

void msm4g_linkedlist_add(LinkedList *list,void *data) {
  LinkedListElement *item;

  item = malloc(sizeof(LinkedListElement));
  item->data = data;
  item->next = NULL;
  item->prev = NULL;

  /* the list is empty */
  if ( list->tail == NULL){
    list->head = item;
    list->tail = item;
  } else {
    item->prev = list->tail;
    list->tail->next = item;
    list->tail = item;
  }
}

void *msm4g_linkedlist_get(LinkedList *list, int index) {
  void *data = NULL;
  LinkedListElement *curr;
  int i;

  curr = list->head;
  i=0;
  while (curr != NULL) {
    if (i == index) {
      data = curr->data;
      break;
    }
    i++;
    curr = curr->next;
  }
  return data;
}

void msm4g_linkedlist_destroy(LinkedList *list) {
  LinkedListElement *current, *next;
  current = list->head;
  while (current != NULL) {
    next=current->next;
    free(current);
    current=next;
  }
  free(list);
  list = NULL;
}

void msm4g_linkedlist_destroyWithData(LinkedList *list) {
  LinkedListElement *curr;
  curr = list->head;
  while (curr != NULL) {
    if (curr->data)
      free(curr->data);
    curr = curr->next;
  }

  msm4g_linkedlist_destroy(list);
}
void msm4g_forwardeuler(double *r,double *v,double *a,double dt,int n,int d,
    double *m,double G) {
  int i;
  msm4g_acceleration(a,r,n,d,m,G);
  for (i=0;i<n*d;i++)
    r[i] += v[i]*dt ;
  for (i=0;i<n*d;i++)
    v[i] += a[i]*dt ;
}

SimulationBox *msm4g_box_new() {
  SimulationBox *box;
  box = malloc(sizeof(SimulationBox));
  msm4g_d3vector_set(&(box->location), 0.0,0.0,0.0);
  msm4g_d3vector_set(&(box->width), 0.0,0.0,0.0);
  return box;
}

SimulationBox *msm4g_box_newCube(double location, double width) {
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

void msm4g_box_update(SimulationBox *box,LinkedList *list,double margin,
    double h,double p) {
  double min[3];
  double max[3];
  double oldwidth,newwidth;
  int i;
  LinkedListElement *curr ;
  Particle *particle;
  for (i=0;i<3;i++) {
    min[i] = DBL_MAX;
    max[i] = DBL_MIN;
  }

  curr = list->head;
  while (curr != NULL) {
    particle = (Particle *) (curr->data);
    for (i=0; i<3; i++) {
      if (particle->r.value[i] >= max[i])
        max[i] = particle->r.value[i];
      if (particle->r.value[i] <= min[i])
        min[i] = particle->r.value[i];
    }
    curr = curr->next;
  }

  /* enclosing the particles tightly */
  for (i=0;i<3;i++) {
    /* If the width in ith dimension is effectively zero,
     * we artificaially add some width (h) so that
     * a grid can be created.
     */
    if (fabs(max[i]-min[i]) < DBL_EPSILON ) {
      box->location.value[i] = min[i]-h*0.50;
      box->width.value[i]    = h;
    } else {
      box->location.value[i] = min[i];
      box->width.value[i] = max[i]-min[i];
    }
  }

  /* enlarging the box <margin> percent. */
  if (margin > 0.0) {
    for (i=0;i<3;i++) {
      oldwidth = box->width.value[i];
      newwidth = oldwidth*(1.0+margin);
      box->location.value[i] -= 0.5*(newwidth-oldwidth);
      box->width.value[i] = newwidth;
    }
  }

  /* make sure the widths are multiple of h */
  for (i=0;i<3;i++) {
    oldwidth = box->width.value[i];
    newwidth = fabs(floor(oldwidth/h)-oldwidth/h) < 0.0001 ? oldwidth :
        (floor(oldwidth/h)+1)*h;
    box->width.value[i] = newwidth;
    box->location.value[i] -= (newwidth-oldwidth)*0.5;
  }

  /* add p/2-1 boxes around the boundary */
  for (i=0;i<3;i++) {
    oldwidth = box->width.value[i];
    newwidth = oldwidth + (p-2)*h ;
    box->width.value[i] = newwidth;
    box->location.value[i] -= h*(p/2-1);
  }
}

void msm4g_box_translate(SimulationBox *box, LinkedList *particles,
    D3Vector delta) {
  int i;
  int n;
  n=3;
  LinkedListElement *curr;
  Particle *particle;

  for (i=0;i<n;i++)
    box->location.value[i] += delta.value[i];
  curr = particles->head;
  while (curr != NULL) {
    particle = (Particle *)curr->data;
    for (i=0;i<n;i++)
      particle->r.value[i] += delta.value[i];
    curr = curr->next;
  }
}

void msm4g_box_translateToOrigin(SimulationBox *box, LinkedList *particles) {
  D3Vector delta;
  delta.value[0] =  -box->location.value[0];
  delta.value[1] =  -box->location.value[1];
  delta.value[2] =  -box->location.value[2];
  msm4g_box_translate(box,particles,delta);
}

void msm4g_box_print(SimulationBox *box) {
  printf("Box Location: %f %f %f\n",
      box->location.value[0],
      box->location.value[1],
      box->location.value[2]);
  printf("Box Width   : %f %f %f\n",
      box->width.value[0],
      box->width.value[1],
      box->width.value[2]);
}

void msm4g_box_destroy(SimulationBox *box) {
  free(box);
}

Bin *msm4g_bin_new(int nx,int ny,int nz) {
  Bin *bin;
  bin = malloc(sizeof(Bin));
  bin->particles = msm4g_linkedlist_new();
  bin->neighbors = msm4g_linkedlist_new();
  bin->nx = nx ;
  bin->ny = ny ;
  bin->nz = nz ;
  bin->cantorindex = msm4g_math_cantor((int [3]){nx,ny,nz}, 3);
  bin->isGhost = false;
  bin->ghostOf = NULL;
  return bin;
}

/** @todo Be carefull, bin indexes can be negative here. */
LinkedList *msm4g_bin_generate(SimulationBox *box,Particle *particles,int n,
    double binwidth,Boolean periodic) {
  LinkedList *binlist;
  Bin *bin;
  binlist = msm4g_linkedlist_new();

  double minrx = DBL_MAX;
  double minry = DBL_MAX;
  double minrz = DBL_MAX;
  double maxrx = DBL_MIN;
  double maxry = DBL_MIN;
  double maxrz = DBL_MIN;

  for (int particleindex = 0 ; particleindex < n ; particleindex++ ) {
    Particle *particle = &(particles[particleindex]);
    double rx = particle->r.value[0];
    double ry = particle->r.value[1];
    double rz = particle->r.value[2];
    int nx = floor(rx/binwidth) + 1;
    int ny = floor(ry/binwidth) + 1;
    int nz = floor(rz/binwidth) + 1;

    if (rx > maxrx ) maxrx = rx;
    if (ry > maxry ) maxry = ry;
    if (rz > maxrz ) maxrz = rz;
    if (rx < minrx ) minrx = rx;
    if (ry < minry ) minry = ry;
    if (rz < minrz ) minrz = rz;

    int binindex = msm4g_math_cantor((int [3]){nx,ny,nz},3);
    bin = msm4g_bin_searchByIndex(binlist,binindex);
    if (bin == NULL) {
      bin = msm4g_bin_new(nx,ny,nz);
      bin->isGhost = false;
      msm4g_linkedlist_add(binlist, bin);
    }
    msm4g_linkedlist_add(bin->particles, particle);
  }

  int minbinindexx = floor(minrx/binwidth) + 1;
  int minbinindexy = floor(minry/binwidth) + 1;
  int minbinindexz = floor(minrz/binwidth) + 1;
  int maxbinindexx = floor(maxrx/binwidth) + 1;
  int maxbinindexy = floor(maxry/binwidth) + 1;
  int maxbinindexz = floor(maxrz/binwidth) + 1;

  msm4g_bin_findneighbors(binlist,minbinindexx,minbinindexy,minbinindexz,
      maxbinindexx,maxbinindexy,maxbinindexz,periodic);

  /* Sanity check: number of particles must be
     equal to sum of particles in the bins */
     int particleCountInBins = 0;
     LinkedListElement *curBin = binlist->head ;
     while (curBin != NULL) {
       Bin *bin = (Bin *)curBin->data;
       if (!bin->isGhost) {
         int count = msm4g_linkedlist_size(bin->particles);
         //msm4g_bin_print(bin);
         particleCountInBins += count;
       }
       curBin = curBin->next;
     }
     if (n != particleCountInBins) {
       fprintf(stderr,"Sanity check in msm4g_bin_generate failed\n");
     }

     //msm4g_bin_printlist(binlist);
     return binlist;
}

void msm4g_bin_findneighbors(LinkedList *binlist,int minbinindexx,
    int minbinindexy, int minbinindexz,int maxbinindexx, int maxbinindexy,
    int maxbinindexz,Boolean periodic) {
  int nbinsx = maxbinindexx - minbinindexx + 1 ;
  int nbinsy = maxbinindexy - minbinindexy + 1 ;
  int nbinsz = maxbinindexz - minbinindexz + 1 ;

  LinkedListElement *curr = binlist->head;
  while (curr != NULL) {
    Bin *bin = (Bin *)curr->data;
    if (bin->isGhost == true) {
      curr = curr->next;
      continue;
    }
    for (int i = -1; i < 2; i++) {
      for (int j = -1; j < 2; j++) {
        for (int k = -1; k < 2; k++) {
          if (i==0 && j==0 && k==0)
            continue;
          int nx = bin->nx + i ;
          int ny = bin->ny + j ;
          int nz = bin->nz + k ;
          int candidateIndex = msm4g_math_cantor((int [3]){nx,ny,nz}, 3);
          Boolean candidateIsGhost = false;
          if (nx < minbinindexx || nx > maxbinindexx ||
              ny < minbinindexy || ny > maxbinindexy ||
              nz < minbinindexz || nz > maxbinindexz)
            candidateIsGhost = true;
          if (!candidateIsGhost) {
            Bin *candidateBin = msm4g_bin_searchByIndex(binlist,candidateIndex);
            if (candidateBin != NULL)
              msm4g_linkedlist_add(bin->neighbors, candidateBin);
          } else if (periodic) {
            Bin *candidateBin = msm4g_bin_searchByIndex(binlist,candidateIndex);
            if (candidateBin != NULL)
              msm4g_linkedlist_add(bin->neighbors, candidateBin);
            else {
              int nxwrapped = nx;
              int nywrapped = ny;
              int nzwrapped = nz;
              if (nxwrapped < minbinindexx) do { nxwrapped += nbinsx ; }
              while ( nxwrapped < minbinindexx) ;
              if (nxwrapped > maxbinindexx) do { nxwrapped -= nbinsx ; }
              while ( nxwrapped > maxbinindexx) ;
              if (nywrapped < minbinindexy) do { nywrapped += nbinsy ; }
              while ( nywrapped < minbinindexy) ;
              if (nywrapped > maxbinindexy) do { nywrapped -= nbinsy ; }
              while ( nywrapped > maxbinindexy) ;
              if (nzwrapped < minbinindexz) do { nzwrapped += nbinsz ; }
              while ( nzwrapped < minbinindexz) ;
              if (nzwrapped > maxbinindexz) do { nzwrapped -= nbinsz ; }
              while ( nzwrapped > maxbinindexz) ;
              int cantorindex = msm4g_math_cantor(
                  (int [3]){nxwrapped,nywrapped,nzwrapped}, 3);
              Bin *ghostOf = msm4g_bin_searchByIndex(binlist, cantorindex);
              if (ghostOf != NULL) {
                Bin *ghostBin = (Bin *)calloc(1,sizeof(Bin));
                ghostBin->particles = ghostOf->particles;
                ghostBin->neighbors = ghostOf->neighbors;
                ghostBin->nx = nx ;
                ghostBin->ny = ny ;
                ghostBin->nz = nz ;
                ghostBin->cantorindex = candidateIndex;
                ghostBin->isGhost = true;
                ghostBin->ghostOf = ghostOf;
                msm4g_linkedlist_add(binlist, ghostBin);
                msm4g_linkedlist_add(bin->neighbors, ghostBin);
              }
            }
          }
        }
      }
    }
    curr = curr->next;
  }
}

Bin *msm4g_bin_searchByIndex(LinkedList *binlist,int cantorindex) {
  LinkedListElement *curr;
  Bin *bin = NULL;
  curr = binlist->head;
  while (curr != NULL) {
    bin = (Bin *)curr->data;
    if (bin->cantorindex == cantorindex) {
      return bin;
    }
    curr = curr->next;
  }

  return NULL;
}

void msm4g_bin_print(Bin *bin) {
  printf("%p [%d,%d,%d] index: %2d ghost:%d\n",
      (void *)bin,bin->nx,bin->ny,bin->nz,
      bin->cantorindex,bin->isGhost);
  //if (!bin->isGhost) {
  LinkedListElement *element = bin->particles->head;

  element = bin->neighbors->head;
  while (element != NULL) {
    Bin *neighbor = (Bin *)element->data;
    printf("  Neighbor %p [%d %d %d  %3d]",element->data,
        neighbor->nx,neighbor->ny,neighbor->nz,neighbor->cantorindex);
    if (neighbor->isGhost)
      printf(" ghost\n");
    else
      printf("\n");
    element = element->next;
  }
}

void msm4g_bin_printlist(LinkedList *binlist) {
  LinkedListElement *curr = binlist->head;
  while (curr != NULL) {
    Bin *bin = (Bin *)curr->data;
    msm4g_bin_print(bin);
    curr = curr->next;
  }
}

void msm4g_bin_destroy(LinkedList *binlist) {
  LinkedListElement *curr;
  Bin *bin;

  curr = binlist->head;
  while (curr != NULL) {
    bin = (Bin *)curr->data;
    if (!bin->isGhost) {
      msm4g_linkedlist_destroy(bin->particles);
      msm4g_linkedlist_destroy(bin->neighbors);
    }
    curr = curr->next;
  }
  msm4g_linkedlist_destroy(binlist);
}

Particle *msm4g_particle_reset(Particle *particle) {
  particle->m = 0.0;
  msm4g_d3vector_set(&(particle->r), 0, 0, 0);
  msm4g_d3vector_set(&(particle->v), 0, 0, 0);
  msm4g_d3vector_set(&(particle->fshort), 0, 0, 0);
  msm4g_d3vector_set(&(particle->flong), 0, 0, 0);
  return particle;
}

Particle **msm4g_particle_rand(int n) {
  Particle **particle;
  int i,j;

  particle = malloc(sizeof(Particle *)*n);

  for (j=0;j<n;j++) {
    particle[j] = msm4g_particle_empty();
    particle[j]->m = (double)rand()/RAND_MAX;
    for (i=0; i<3; i++) {
      particle[j]->r.value[i] = (double)rand()/RAND_MAX;
      particle[j]->v.value[i] = (double)rand()/RAND_MAX;
      particle[j]->fshort.value[i] = (double)rand()/RAND_MAX;
      particle[j]->flong.value[i] = (double)rand()/RAND_MAX;
    }
  }
  return particle;
}

Particle *msm4g_particle_read(const char *filename,int *numberofparticles) {
  Particle *particles;
  FILE *fp;
  int n;
  double mass;
  double r[3];
  double v[3];

  fp = fopen(filename,"r");
  if (fp == NULL) {
    return NULL;
  }
  fscanf(fp,"%d",&n);
  *numberofparticles = n;
  particles = (Particle *)calloc(n,sizeof(Particle));
  for (int iparticle = 0 ; iparticle < n ; iparticle++) {
    int ret = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
        &mass,&r[0],&r[1],&r[2],&v[0],&v[1],&v[2]);
    if (ret == 7){
      particles[iparticle].m = mass;
      particles[iparticle].index = iparticle;
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

Particle *msm4g_particle_empty() {
  static int index = 0;
  int i;
  Particle *particle;

  particle = malloc(sizeof(Particle));
  particle->index = index;
  particle->m = 0.0;
  for (i=0; i<3; i++) {
    particle->r.value[i] = 0.0;
    particle->v.value[i] = 0.0;
    particle->fshort.value[i] = 0.0;
    particle->flong.value[i] = 0.0;
    particle->acc_short[i] = 0.0;

  }

  particle->potential_short_real = 0.0;
#ifdef DEBUG
  particle->potential_short_real_true = 0.0;
#endif
  index++;
  return particle;
}

Particle *msm4g_particle_new(double mass,double x,double y,double z) {
  Particle *particle;
  particle = msm4g_particle_empty();
  particle->m = mass;
  particle->r.value[0] = x;
  particle->r.value[1] = y;
  particle->r.value[2] = z;
  return particle;
}

void msm4g_particle_print(Particle *particle) {
  printf("r:%8.2E %8.2E %8.2E ",particle->r.value[0],
      particle->r.value[1],particle->r.value[2]);
}

void msm4g_particle_printlist(LinkedList *particlelist) {
  Particle *particle;
  LinkedListElement *curr;
  int i=0;
  curr = particlelist->head;
  while (curr != NULL) {
    particle = curr->data;
    printf("Particle index: %d\n",i);
    msm4g_particle_print(particle);

    i++;
    curr=curr->next;
  }
}

void msm4g_particle_destroy(Particle *particle) {
  free(particle);
  particle=NULL;
}

void msm4g_particle_destroyarray(Particle **particlearray,int length) {
  int i;
  for (i=0; i<length; i++) {
    free(particlearray[i]);
    particlearray[i] = NULL;
  }
  free(particlearray);
  particlearray = NULL;
}

void msm4g_acceleration(double *a,double *r,int n,int d,double *m,double G) {
  int i,j,k;
  double rij2,rij3,acc;
  msm4g_zeros(a,n*d);
  for (i=0;i<n-1;i++) {
    for (j=i+1;j<n;j++) {
      rij2=0.0;
      for (k=0;k<d;k++)
        rij2 +=(r[d*i+k]-r[d*j+k])*(r[d*i+k]-r[d*j+k]);
      rij3=sqrt(rij2)*rij2;
      for (k=0;k<d;k++) {
        acc=G*(r[d*i+k]-r[d*j+k])/rij3;
        a[d*j+k] += m[i]*acc;
        a[d*i+k] -= m[j]*acc;
      }
    }
  }
}


void msm4g_energy(double *pot,double *kin,double *tot,double *r,double *v,
    int n,int d,double *m,double G) {
  int i,j,k;
  double rij2,v2;
  *pot=0; *kin=0; *tot=0;
  for (i=0;i<n-1;i++) {
    for (j=i+1;j<n;j++) {
      rij2=0;
      for (k=0;k<d;k++) {
        rij2 += pow(r[d*i+k]-r[d*j+k],2);
      }
      *pot -= m[i]*m[j]/sqrt(rij2);
    }
  }
  *pot *= G;

  for (i=0;i<n;i++) {
    v2 = 0;
    for (k=0;k<d;k++) {
      v2 += v[d*i+k]*v[d*i+k];
    }
    *kin += m[i]*v2;
  }
  *kin *= 0.5;
  *tot = (*pot) + (*kin);
}



void msm4g_print_energy(double pot0,double kin0,double tot0,double pot,
    double kin,double tot) {
  printf("[pot: %10.3E] [kin: %10.3E] [tot: %10.3E] [err: %10.3E]\n",
      pot,kin,tot,(tot-tot0)/tot0);
}

void msm4g_leapfrog(double *r,double *v,double *a,double *a1,double dt,
    int n,int d,double *m,double G) {
  int i;
  msm4g_acceleration(a,r,n,d,m,G);
  for (i=0;i<n*d;i++)
    r[i] += v[i]*dt + a[i]*0.5*dt*dt;
  msm4g_acceleration(a1,r,n,d,m,G);
  for (i=0;i<n*d;i++)
    v[i] += (a[i]+a1[i])*0.5*dt;
}

void msm4g_zeros(double *x,int n) {
  int i;
  for (i=0;i<n;i++)
    x[i]=0.0;
}

int msm4g_math_cantor(int *x,int n) {
  int i,j,h;
  int hsum = 0;
  int jprod = 1;
  int isum=0;
  for (h=1; h<=n; h++) {
    jprod = 1;
    for (j=0; j<= h-1; j++) {
      isum = 0;
      for (i=1; i<=h; i++) {
        isum += x[i-1];
      }
      jprod *= isum + j;
    }
    hsum += jprod/msm4g_math_factorial(h);
  }
  return hsum;
}

int msm4g_math_factorial(int n) {
  int factorial = 1;
  int i;
  if (n==0) return 1;
  for (i=1; i<= n; i++) {
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
    double relerr = msm4g_util_diffnorm(g, gold, p / 2)/msm4g_util_norm(g,p/2);
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

  for (int i = 0; i < p - 1; i++) {
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

int msm4g_util_face_enumerate(int n,SimulationParameters *sp) {
  int face_len = 0 ;
  int counter = 0;
  if (n == 0) {
    sp->face_i[counter] = 0;
    sp->face_j[counter] = 0;
    sp->face_k[counter] = 0;
    counter++;
  } else {
    for (int j = -n ; j <= n ; j++) {
      for (int k = -n ; k <= n ; k++) {
        sp->face_i[counter] = -n ;
        sp->face_j[counter] =  j ;
        sp->face_k[counter] =  k  ;
        counter++;
        sp->face_i[counter] = +n ;
        sp->face_j[counter] =  j ;
        sp->face_k[counter] =  k  ;
        counter++;
      }
    }
    for (int i = -n+1 ; i < n ; i++) {
      for (int k = -n ; k <= n ; k++) {
        sp->face_i[counter] =  i ;
        sp->face_j[counter] = -n ;
        sp->face_k[counter] =  k ;
        counter++;
        sp->face_i[counter] =  i ;
        sp->face_j[counter] =  n ;
        sp->face_k[counter] =  k ;
        counter++;
      }
    }
    for (int i = -n+1 ; i < n ; i++) {
      for (int j = -n+1 ; j < n ; j++) {
        sp->face_i[counter] =  i ;
        sp->face_j[counter] =  j ;
        sp->face_k[counter] = -n ;
        counter++;
        sp->face_i[counter] =  i ;
        sp->face_j[counter] =  j ;
        sp->face_k[counter] =  n ;
        counter++;
      }
    }
  }
  face_len = counter;
  if (counter >= FACE_MAXLEN) {
    fprintf(stderr,"face arrays are not large enough\n");
    face_len = FACE_MAXLEN - 1 ;
  }
  return face_len ;
}

double msm4g_util_choose_beta(double aL,double tol,double h0) {
  double const_ = tol*aL/h0;
  double beta = 1;
  double res = erfc(beta) - const_;
  double oldres = 2.*res;
  while (fabs(res) < fabs(oldres)){
    double dres = -2./sqrt(MYPI)*exp(-beta*beta);
    beta -= res/dres;
    oldres = res;
    res = erfc(beta) - const_;
  }
  beta = beta/aL;
  return beta;
}

double msm4g_util_choose_kmax(double beta,double tol,double h0) {  
  double const_ = sqrt(MYPI)*tol/(2.*beta*h0);
  double kmax = 1.;
  double res = erfc(kmax) - const_;
  double oldres = 2.*res;
  while (fabs(res) < fabs(oldres)){
    double dres = -2./sqrt(MYPI)*exp(-kmax*kmax);
    kmax -= res/dres;
    oldres = res;
    res = erfc(kmax) - const_;}
  kmax = beta*kmax/MYPI;
  return kmax;
}

double msm4g_util_calculate_c(int k, double M, int nu, int mu, double *wprime) {
  double c = wprime[0];
  // Use the fact that sin components cancel
  for (int m = 1; m <= mu + nu/2; m++) {
    c += 2.0 * wprime[m]*cos(2.0 * MYPI * k * m / (double)M) ;
  }
  return c;
}

double msm4g_util_nchoosek(int n,int k) {
  double out = 1.0;
  if (k < 0 || k > n) {
    out = 0.0;
  } else if (k == 0 ) {
    out = 1.0;
  } else {
    for (int i = 0 ; i < k ;i++)
      out *= (n-i) / (double)(i+1) ;
  }
  return out;
}

double msm4g_util_jn(int nu,int k) {
  return pow(2,1-nu) * msm4g_util_nchoosek(nu, nu/2 + abs(k));
}
