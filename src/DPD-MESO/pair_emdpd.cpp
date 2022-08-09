// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Zhen Li (Clemson University)
   Email: zli7@clemson.edu
------------------------------------------------------------------------- */

#include "pair_emdpd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define EPSILON 1.0e-10

static const char cite_pair_emdpd[] =
  "pair emdpd command:\n\n"
  "@Article{ZLi2014_JCP,\n"
  " author = {Li, Z. and Tang, Y.-H. and Lei, H. and Caswell, B. and Karniadakis, G.E.},\n"
  " title = {Energy-conserving dissipative particle dynamics with temperature-dependent properties},\n"
  " journal = {Journal of Computational Physics},\n"
  " year =    {2014},\n"
  " volume =  {265},\n"
  " pages =   {113--127}\n"
  "}\n\n"
  "@Article{ZLi2015_CC,\n"
  " author = {Li, Z. and Tang, Y.-H. and Li, X. and Karniadakis, G.E.},\n"
  " title = {Mesoscale modeling of phase transition dynamics of thermoresponsive polymers},\n"
  " journal = {Chemical Communications},\n"
  " year =    {2015},\n"
  " volume =  {51},\n"
  " pages =   {11038--11040}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairEMDPD::PairEMDPD(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_emdpd);
  writedata = 1;
  random = nullptr;
  randomT = nullptr;
}

/* ---------------------------------------------------------------------- */

PairEMDPD::~PairEMDPD()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_r);
    memory->destroy(cutT);

    memory->destroy(A_att);
    memory->destroy(B_rep);
    memory->destroy(gamma);
    memory->destroy(power);
    memory->destroy(kappa);
    memory->destroy(powerT);
  }
  if (power_flag) memory->destroy(sc);
  if (kappa_flag) memory->destroy(kc);

  if (random) delete random;
  if (randomT) delete randomT;
}

/* ---------------------------------------------------------------------- */

void PairEMDPD::compute(int eflag, int vflag)
{
  double evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *T = atom->emdpd_temp;
  double *Q = atom->emdpd_flux;
  double *cv = atom->emdpd_cv;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  double *rho = atom->rho;
  double *phi = atom->phi;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);
  double kboltz = 1.0;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    double vxtmp = v[i][0];
    double vytmp = v[i][1];
    double vztmp = v[i][2];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    double rhoi = rho[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      double factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        double r = sqrt(rsq);
        if (r < EPSILON) continue;
        double rinv = 1.0/r;
        double delvx = vxtmp - v[j][0];
        double delvy = vytmp - v[j][1];
        double delvz = vztmp - v[j][2];
        double dot = delx*delvx + dely*delvy + delz*delvz;
        double vijeij = dot*rinv;
        double randnum = random->gaussian();

        double T_ij=0.5*(T[i]+T[j]);
        double T_pow[4];
        T_pow[0] = T_ij - 1.0;
        T_pow[1] = T_pow[0]*T_pow[0];
        T_pow[2] = T_pow[0]*T_pow[1];
        T_pow[3] = T_pow[0]*T_pow[2];

        double power_d = power[itype][jtype];
        if (power_flag) {
          double factor = 1.0;
          for (int k = 0; k < 4; k++)
            factor += sc[itype][jtype][k]*T_pow[k];
          power_d *= factor;
        }

        power_d = MAX(0.01,power_d);
        double wc = 1.0 - r/cut[itype][jtype];
        wc = MAX(0.0,MIN(1.0,wc));
        double wc_r = 1.0 - r/cut_r[itype][jtype];
        wc_r = MAX(0.0,wc_r);
        double wr = pow(wc, 0.5*power_d);

        double GammaIJ = gamma[itype][jtype];
        double SigmaIJ = 4.0*GammaIJ*kboltz*T[i]*T[j]/(T[i]+T[j]);
        SigmaIJ = sqrt(SigmaIJ);
        double rhoj = rho[j];

        double ratio = 1.0;
        if(itype != flag_wall && jtype == flag_wall)
          ratio = effective_factor(phi[i],cut[itype][jtype]);
        else if(jtype != flag_wall && itype == flag_wall) 
          ratio = effective_factor(phi[j],cut[itype][jtype]);
   
        SigmaIJ *= sqrt(ratio);
        GammaIJ *= ratio;

        double fpair =  A_att[itype][jtype]*T_ij*wc + B_rep[itype][jtype]*(rhoi+rhoj)*wc_r;
        fpair -= GammaIJ *wr*wr *dot*rinv;
        fpair += SigmaIJ * wr *randnum * dtinvsqrt;
        fpair *= factor_dpd*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        // heat transfer
        double dQc,dQd,dQr;
        if (r < cutT[itype][jtype]) {
          double wrT = 1.0 - r/cutT[itype][jtype];
          wrT = MAX(0.0,MIN(1.0,wrT));
          wrT = pow(wrT, 0.5*powerT[itype][jtype]);
          double randnumT = randomT->gaussian();
          randnumT = MAX(-5.0,MIN(randnum,5.0));

          double kappaT = kappa[itype][jtype];
          if (kappa_flag) {
            double factor = 1.0;
            for (int k = 0; k < 4; k++)
              factor += kc[itype][jtype][k]*T_pow[k];
            kappaT *= factor;
          }

          double kij = cv[i]*cv[j]*kappaT * T_ij*T_ij;
          double alphaij = sqrt(2.0*kboltz*kij);

          dQc  = kij * wrT*wrT * ( T[j] - T[i] )/(T[i]*T[j]);
          dQd  = wr*wr*( GammaIJ * vijeij*vijeij - SigmaIJ*SigmaIJ/mass[itype] ) - SigmaIJ * wr *vijeij *randnum;
          dQd /= (cv[i]+cv[j]);
          dQr  = alphaij * wrT * dtinvsqrt * randnumT;
          Q[i] += (dQc + dQd + dQr );
        }
        //-----------------------------------------------------------

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
          Q[j]    -= ( dQc - dQd + dQr );
        }

        if (eflag) {
          evdwl = 0.5*A_att[itype][jtype]*T_ij*cut[itype][jtype] * wc*wc;
          evdwl += 0.5*B_rep[itype][jtype]*cut_r[itype][jtype]*(rhoi+rhoj)*wc_r*wc_r;
          evdwl *= factor_dpd;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEMDPD::allocate()
{
  int i,j;
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (i = 1; i <= n; i++)
  for (j = i; j <= n; j++)
    setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_r,n+1,n+1,"pair:cut_r");
  memory->create(cutT,n+1,n+1,"pair:cutT");
  memory->create(A_att,n+1,n+1,"pair:A_att");
  memory->create(B_rep,n+1,n+1,"pair:B_rep");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(power,n+1,n+1,"pair:power");
  memory->create(kappa,n+1,n+1,"pair:kappa");
  memory->create(powerT,n+1,n+1,"pair:powerT");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEMDPD::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  seed = utils::inumeric(FLERR,arg[1],false,lmp);
  flag_wall = utils::inumeric(FLERR,arg[2],false,lmp);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0)
    error->all(FLERR,"Invalid random number seed");

  delete random;
  random = new RanMars(lmp,(seed + comm->me) % 900000000);
  randomT = new RanMars(lmp,(2*seed + comm->me) % 900000000);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
    for (j = i+1; j <= atom->ntypes; j++)
    if (setflag[i][j])
      cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEMDPD::coeff(int narg, char **arg)
{
  if (narg < 11)
    error->all(FLERR,"Incorrect args for pair emdpd coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double A_one = utils::numeric(FLERR,arg[2],false,lmp);
  double B_one = utils::numeric(FLERR,arg[3],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[4],false,lmp);
  double power_one = utils::numeric(FLERR,arg[5],false,lmp);
  double cut_one   = utils::numeric(FLERR,arg[6],false,lmp);
  double cut_two   = utils::numeric(FLERR,arg[7],false,lmp);
  double kappa_one = utils::numeric(FLERR,arg[8],false,lmp);
  double powerT_one= utils::numeric(FLERR,arg[9],false,lmp);
  double cutT_one  = utils::numeric(FLERR,arg[10],false,lmp);

  if (cut_one < cut_two) error->all(FLERR,"Incorrect args for pair coefficients\n cutA should be larger than cutB.");

  int iarg = 11;
  power_flag = kappa_flag = 0;
  double sc_one[4], kc_one[4];
  int n = atom->ntypes;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"power") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal pair emdpd coefficients");
      for (int i = 0; i < 4; i++)
        sc_one[i] = utils::numeric(FLERR,arg[iarg+i+1],false,lmp);
      iarg += 5;
      power_flag = 1;
      memory->create(sc,n+1,n+1,4,"pair:sc");
    } else if (strcmp(arg[iarg],"kappa") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal pair emdpd coefficients");
      for (int i = 0; i < 4; i++)
        kc_one[i] = utils::numeric(FLERR,arg[iarg+i+1],false,lmp);
      iarg += 5;
      kappa_flag = 1;
      memory->create(kc,n+1,n+1,4,"pair:kc");
    } else error->all(FLERR,"Illegal pair emdpd coefficients");
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++)
  for (int j = MAX(jlo,i); j <= jhi; j++) {
    A_att[i][j] = A_one;
    B_rep[i][j] = B_one;
    gamma[i][j] = gamma_one;
    power[i][j] = power_one;
    cut[i][j]   = cut_one;
    cut_r[i][j] = cut_two;
    kappa[i][j] = kappa_one;
    powerT[i][j]= powerT_one;
    cutT[i][j]  = cutT_one;

    if (power_flag)
    for (int k = 0; k < 4; k++)
      sc[i][j][k] = sc_one[k];

    if (kappa_flag)
    for (int k = 0; k < 4; k++)
      kc[i][j][k] = kc_one[k];

    setflag[i][j] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEMDPD::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair emdpd requires ghost atoms store velocity");

  if (!atom->rho_flag)
    error->all(FLERR,"Pair style mdpd requires atom attribute rho");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR, "Pair tdpd needs newton pair on for momentum conservation");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEMDPD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cut[j][i]   = cut[i][j];
  cut_r[j][i] = cut_r[i][j];
  cutT[j][i]  = cutT[i][j];
  A_att[j][i] = A_att[i][j];
  B_rep[j][i] = B_rep[i][j];
  gamma[j][i] = gamma[i][j];
  power[j][i] = power[i][j];
  kappa[j][i] = kappa[i][j];
  powerT[j][i]= powerT[i][j];

  if (power_flag)
  for (int k = 0; k < 4; k++)
    sc[j][i][k] = sc[i][j][k];

  if (kappa_flag)
  for (int k = 0; k < 4; k++)
    kc[j][i][k] = kc[i][j][k];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEMDPD::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
  for (int j = i; j <= atom->ntypes; j++) {
    fwrite(&setflag[i][j],sizeof(int),1,fp);
    if (setflag[i][j]) {
      fwrite(&A_att[i][j],sizeof(double),1,fp);
      fwrite(&B_rep[i][j],sizeof(double),1,fp);
      fwrite(&gamma[i][j],sizeof(double),1,fp);
      fwrite(&power[i][j],sizeof(double),1,fp);
      fwrite(&cut[i][j],sizeof(double),1,fp);
      fwrite(&cut_r[i][j],sizeof(double),1,fp);
      fwrite(&kappa[i][j],sizeof(double),1,fp);
      fwrite(&powerT[i][j],sizeof(double),1,fp);
      fwrite(&cutT[i][j],sizeof(double),1,fp);
      if (power_flag)
      for (int k = 0; k < 4; k++)
        fwrite(&sc[i][j][k],sizeof(double),1,fp);

      if (kappa_flag)
      for (int k = 0; k < 4; k++)
        fwrite(&kc[i][j][k],sizeof(double),1,fp);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEMDPD::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&A_att[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&B_rep[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&power[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_r[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&kappa[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&powerT[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cutT[i][j],sizeof(double),1,fp,nullptr,error);
          if (power_flag)
          for (int k = 0; k < 4; k++)
            utils::sfread(FLERR,&sc[i][j][k],sizeof(double),1,fp,nullptr,error);

          if (kappa_flag)
          for (int k = 0; k < 4; k++)
            utils::sfread(FLERR,&kc[i][j][k],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&A_att[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&B_rep[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&power[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_r[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&kappa[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&powerT[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutT[i][j],1,MPI_DOUBLE,0,world);
        if (power_flag)
        for (int k = 0; k < 4; k++)
          MPI_Bcast(&sc[i][j][k],1,MPI_DOUBLE,0,world);

        if (kappa_flag)
        for (int k = 0; k < 4; k++)
          MPI_Bcast(&kc[i][j][k],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEMDPD::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEMDPD::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
  if (randomT) delete randomT;
  randomT = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

double PairEMDPD::single(int i, int j, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_dpd, double &fforce)
{
  double r,rinv,wc,phi;
  double *T = atom->emdpd_temp;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }
  double T_ij = 0.5*(T[i]+T[j]);
  rinv = 1.0/r;
  wc = 1.0 - r/cut[itype][jtype];
  fforce = A_att[itype][jtype]*T_ij*wc*factor_dpd*rinv;

  phi = 0.5*A_att[itype][jtype]*T_ij*cut[itype][jtype]*wc*wc;
  return factor_dpd*phi;
}

double PairEMDPD::effective_factor(double phi, double rcw)
{
  double h, h_inv, ratio;
  double rcw_inv = 1.0/rcw;
  phi = phi < 0.5 ? phi : 0.5; 
  
  h= 1 - pow(2.088*phi*phi*phi + 1.478*phi, 0.25);
  h = h > 0.025 ? h : 0.025;

  h *= rcw;
  h_inv = 1.0/h;

  ratio = 1.0 + 0.187*(rcw*h_inv - 1.0) - 0.093*(1.0-h*rcw_inv)*(1.0-h*rcw_inv)*(1.0-h*rcw_inv);

  return ratio;
}
