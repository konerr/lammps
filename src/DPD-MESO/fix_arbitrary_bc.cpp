/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This is a fix for implementation of Arbitrary BC combined with 
   the density computation.

   Contributing author: Zhen Li (Clemson University)
   Email: zli7@clemson.edu
------------------------------------------------------------------------- */

#include <mpi.h>
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_arbitrary_bc.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair_emdpd.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "citeme.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#define const_pi 3.141592653589793
#define EPSILON 1.0e-6

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_arbitrarybc[] =
  "fix arbitrary/bc command:\n\n"
  "@Article{ZLi2017_JCP,\n"
  " author = {Li, Z. and Bian, X. and Tang, Y.-H. and Karniadakis, G.E.},\n"
  " title = {A dissipative particle dynamics method for arbitrarily complex geometries},\n"
  " journal = {Journal of Computational Physics},\n"
  " year = {2018},\n"
  " volume = {355},\n"
  " pages = {534-547}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixArbitraryBc::FixArbitraryBc(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  emdpd(NULL)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_arbitrarybc);

  if (strcmp(style,"abitrary/bc") != 0 && narg < 8)
    error->all(FLERR,"Illegal fix abitrary/bc command: flag_wall rho_wall cut_phi phi_c cut_rho");

  /* flag_wall = force->inumeric(FLERR,arg[3]); */
  /* rho_wall = force->numeric(FLERR,arg[4]); */
  /* cut_phi = force->numeric(FLERR,arg[5]); */
  /* phi_c = force->numeric(FLERR,arg[6]); */
  /* cut_rho = force->numeric(FLERR,arg[7]); */

  flag_wall = utils::inumeric(FLERR,arg[3],false,lmp);
  rho_wall  = utils::numeric(FLERR,arg[4],false,lmp);
  cut_phi   = utils::numeric(FLERR,arg[5],false,lmp);
  phi_c     = utils::numeric(FLERR,arg[6],false,lmp);
  cut_rho   = utils::numeric(FLERR,arg[7],false,lmp);
  
  newton_pair = force->newton_pair;
  if (newton_pair) comm_reverse = 1;
   
  double rc_3 = cut_rho*cut_rho*cut_rho;
  double rc_7 = rc_3*rc_3*cut_rho;
  rho_factor = 105.0 / (16.0* const_pi *rc_7);

  rc_3 = cut_phi*cut_phi*cut_phi;
  rc_7 = rc_3*rc_3*cut_phi;
  phi_factor = 105.0 / (16.0* const_pi * rho_wall *rc_7);
  dw_factor = -315.0 / (4.0 * const_pi * rho_wall *rc_7);

  emdpd = NULL;
  comm_reverse = 5;
}

/* ---------------------------------------------------------------------- */

FixArbitraryBc::~FixArbitraryBc()
{
}

/* ---------------------------------------------------------------------- */

int FixArbitraryBc::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixArbitraryBc::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  emdpd = (PairEMDPD *) force->pair_match("emdpd",1);
  if (emdpd == NULL)
    error->all(FLERR,"Must use pair_style emdpd with fix arbitrary/bc");
}

void FixArbitraryBc::post_integrate()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double r,rsq,h,r_inv,h_inv,wf;
  double vel_half[3], dot, tmp, dtfm, norm_inv;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho = atom->rho;
  double *phi= atom->phi;
  double **nw = atom->nw;
  double **vest = atom->vest;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  NeighList *list = emdpd->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int np = nlocal;
  if(newton_pair) np = nall;

  for( i = 0; i < np; i++) {
    rho[i] = 0.0;
    phi[i] = 0.0; 
    for( j = 0; j < 3; j++) {
      nw[i][j] = 0.0;
    }
  }    

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      r = sqrt(rsq);
      r = MAX(r,EPSILON);

      h = cut_rho;
      if (r < h) {
        wf = (h + 3.0*r)*(h - r)*(h - r)*(h - r)*rho_factor;
        rho[i] += mass[jtype]*wf;
        if (newton_pair || j < nlocal)
          rho[j] += mass[itype]*wf;
      }

      h = cut_phi;
      if ( (itype==flag_wall || jtype==flag_wall) && r < h ) {
        r_inv = 1.0/r;
        if(itype == flag_wall && jtype != flag_wall) {
          phi[j] += (h + 3.0*r)*(h - r)*(h - r)*(h - r)*phi_factor;
          tmp = r*(h-r)*(h-r)*dw_factor;
          nw[j][0] += tmp * delx *r_inv;
          nw[j][1] += tmp * dely *r_inv;
          nw[j][2] += tmp * delz *r_inv;
        }
        else if(itype != flag_wall && jtype == flag_wall) {
          phi[i] += (h + 3.0*r)*(h - r)*(h - r)*(h - r)*phi_factor;
          tmp = r*(h-r)*(h-r)*dw_factor;
          nw[i][0] -= tmp * delx *r_inv;
          nw[i][1] -= tmp * dely *r_inv;
          nw[i][2] -= tmp * delz *r_inv;    
        }
      }
    }
  }
  /* if (newton_pair) comm->reverse_comm_fix_variable(this); */ 
  if (newton_pair) comm->reverse_comm(this); 
 
  for (i = 0; i < nlocal; i++)
  if ( (mask[i] & groupbit) && phi[i] > phi_c) {
     dtfm = dtf / mass[type[i]];

     norm_inv = 1.0/sqrt(nw[i][0]*nw[i][0] + nw[i][1]*nw[i][1] + nw[i][2]*nw[i][2]);
     dot = 0.0;
     for(int k = 0; k < 3; k++) {
       nw[i][k] *= norm_inv;
       vel_half[k] = vest[i][k] + dtfm * f[i][k];
       x[i][k] -= dtv * vel_half[k];
       dot += vel_half[k]*nw[i][k];
     }
     
     tmp = MAX(0.0,dot);
     for(int k = 0; k < 3; k++) {     
       v[i][k] = -vel_half[k] + 2.0*tmp*nw[i][k];
       x[i][k] += dtv * v[i][k];
     }
  }
}

int FixArbitraryBc::pack_reverse_comm(int n, int first, double *buf) {
  int i, m, last;
  double *rho = atom->rho;
  double *phi = atom->phi;
  double **nw = atom->nw;  

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = rho[i];
    buf[m++] = phi[i];
    buf[m++] = nw[i][0];
    buf[m++] = nw[i][1];
    buf[m++] = nw[i][2];
  }
  return m;
}

void FixArbitraryBc::unpack_reverse_comm(int n, int *list, double *buf) {
  int i, j, m;
  double *rho = atom->rho;
  double *phi = atom->phi;
  double **nw = atom->nw;  

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
    phi[j] += buf[m++];
    nw[j][0] += buf[m++];
    nw[j][1] += buf[m++];
    nw[j][2] += buf[m++];
  }
}

