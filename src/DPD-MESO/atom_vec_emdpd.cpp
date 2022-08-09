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

#include "atom_vec_emdpd.h"

#include "atom.h"
#include "error.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecEMDPD::AtomVecEMDPD(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->emdpd_flag = 1;
  atom->vest_flag = 1;
  atom->rho_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"emdpd_cv", "emdpd_temp", "emdpd_flux", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_copy = {"emdpd_cv", "emdpd_temp", "emdpd_flux", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_comm = {"emdpd_temp", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_comm_vel = {"emdpd_temp", "vest", "vest_temp"};
  fields_reverse = {"emdpd_flux"};
  fields_border = {"emdpd_cv", "emdpd_temp", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_border_vel = {"emdpd_cv", "emdpd_temp", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_exchange = {"emdpd_cv", "emdpd_temp", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_restart = {"emdpd_cv", "emdpd_temp", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_create = {"emdpd_cv", "emdpd_temp", "emdpd_flux", "vest", "vest_temp", "rho", "phi", "nw"};
  fields_data_atom = {"id", "type", "emdpd_temp", "emdpd_cv", "x", "rho", "phi", "nw"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecEMDPD::init()
{
  AtomVec::init();

  if (strcmp(update->unit_style, "lj") != 0) error->all(FLERR, "Atom style emdpd requires lj units");
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecEMDPD::grow_pointers()
{
  emdpd_cv = atom->emdpd_cv;
  emdpd_temp = atom->emdpd_temp;
  emdpd_flux = atom->emdpd_flux;
  vest = atom->vest;
  vest_temp = atom->vest_temp;
  rho = atom->rho;
  phi = atom->phi;
  nw = atom->nw;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecEMDPD::force_clear(int n, size_t nbytes)
{
  memset(&emdpd_flux[n], 0, nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecEMDPD::create_atom_post(int ilocal)
{
  emdpd_temp[ilocal] = 1.0;
  emdpd_cv[ilocal] = 1.0e5;
  vest_temp[ilocal] = emdpd_temp[ilocal];
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecEMDPD::data_atom_post(int ilocal)
{
  emdpd_flux[ilocal] = 0.0;
  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;
  vest_temp[ilocal] = emdpd_temp[ilocal];
}
