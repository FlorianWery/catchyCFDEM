/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(read_dump_to_fix,ReadDumpToFix)

#else

#ifndef LMP_READ_DUMP_TO_FIX_H
#define LMP_READ_DUMP_TO_FIX_H

#include <stdio.h>
#include "pointers.h"

namespace LAMMPS_NS {

class ReadDumpToFix : protected Pointers {
 friend class Rerun;

 public:
  ReadDumpToFix(class LAMMPS *);
  ~ReadDumpToFix();
  void command(int, char **);

  void store_files(int, char **);
  void setup_reader(int, char **);
  bigint seek(bigint, int);
  void header(int);
  bigint next(bigint, bigint, int, int);
  void atoms();
  int fields_and_keywords(int, char **);

private:

  void file_search(char *infile);

  int me,nprocs;
  FILE *fp;

  int dimension;
  int triclinic;

  int nfile;               // # of dump files to process
  char **files;            // list of file names
  int currentfile;         // currently open file

  int boxflag;             // overwrite simulation with dump file box params
  int replaceflag,addflag; // flags for processing dump snapshot atoms
  int trimflag,purgeflag;
  int scaleflag;           // user 0/1 if dump file coords are unscaled/scaled
  int wrapflag;            // user 0/1 if dump file coords are unwrapped/wrapped
  char *readerstyle;       // style of dump files to read

  int nfield;              // # of fields to extract from dump file
  int *fieldtype;          // type of each field = X,VY,IZ,etc
  char **fieldlabel;       // user specified label for field
  double **fields;         // per-atom field values

  int scaled;              // 0/1 if dump file coords are unscaled/scaled
  int wrapped;             // 0/1 if dump file coords are unwrapped/wrapped

  double box[3][3];         // dump file box parameters
  double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz;  // dump snapshot box params
  double xprd,yprd,zprd;

  bigint nsnapatoms;        // # of atoms in dump file shapshot

  int npurge,nreplace,ntrim,nadd;     // stats on processed atoms
  int addproc;                        // proc that should add next atom
  int yindex,zindex;                  // field index for Y,Z coords

  int *uflag;               // set to 1 if snapshot atom matches owned atom
  int *ucflag,*ucflag_all;  // set to 1 if snapshot chunk atom was processed

  class Reader *reader;           // class that reads dump file

  void process_atoms(int);
  void delete_atoms();

  double xfield(int, int);
  double yfield(int, int);
  double zfield(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Read_dump command before simulation box is defined

The read_dump command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump file does not contain requested snapshot

Self-explanatory.

E: Invalid dump reader style

Self-explanatory.

E: No box information in dump. You have to use 'box no'

Self-explanatory.

E: Read_dump triclinic status does not match simulation

Both the dump snapshot and the current LAMMPS simulation must
be using either an orthogonal or triclinic box.

E: Read_dump field not found in dump file

Self-explanatory.

E: Read_dump x,y,z fields do not have consistent scaling

Self-explanatory.

E: All read_dump x,y,z fields must be specified for scaled, triclinic coords

For triclinic boxes and scaled coordinates you must specify all 3 of
the x,y,z fields, else LAMMPS cannot reconstruct the unscaled
coordinates.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: Duplicate fields in read_dump command

Self-explanatory.

E: If read_dump purges it cannot replace or trim

These operations are not compatible.  See the read_dump doc
page for details.

*/
