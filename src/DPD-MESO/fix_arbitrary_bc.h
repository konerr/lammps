/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(arbitrary/bc,FixArbitraryBc)

#else

#ifndef LMP_FIX_ARBITRARY_BC_H
#define LMP_FIX_ARBITRARY_BC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixArbitraryBc : public Fix {
 public:
  FixArbitraryBc(class LAMMPS *, int, char **);
  ~FixArbitraryBc();
  int setmask();
  virtual void init();
  virtual void post_integrate();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 protected:
  double phi_c, rho_wall, cut_phi, cut_rho;
  int flag_wall;
  int newton_pair;
  double rho_factor, phi_factor, dw_factor;
  double dtv,dtf;
 
  class PairEMDPD *emdpd;
};

}

#endif
#endif
