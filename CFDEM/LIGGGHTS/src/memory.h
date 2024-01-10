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

#ifndef LMP_MEMORY_H
#define LMP_MEMORY_H

#include "lmptype.h"
#include "pointers.h"

#include <iostream>

namespace LAMMPS_NS {

class Memory : protected Pointers {
 public:
  Memory(class LAMMPS *);

  void *smalloc(bigint n, const char *);
  void *srealloc(void *, bigint n, const char *);
  void sfree(void *);
  void fail(const char *);

/* ----------------------------------------------------------------------
   create/grow/destroy vecs and multidim arrays with contiguous memory blocks
   only use with primitive data types, e.g. 1d vec of ints, 2d array of doubles
   cannot use with pointers, e.g. 1d vec of int*, due to mismatched destroy
   avoid use with non-primitive data types to avoid code bloat
   for these other cases, use smalloc/srealloc/sfree directly
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * n;
      array = (TYPE *) smalloc(nbytes,name);
      return array;
    }

  template <typename TYPE>
    TYPE **create(TYPE **&array, int n, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   grow or shrink 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *grow(TYPE *&array, int n, const char *name)
    {
      if (array == NULL) return create(array,n,name);

      bigint nbytes = ((bigint) sizeof(TYPE)) * n;
      array = (TYPE *) srealloc(array,nbytes,name);
      return array;
    }

  template <typename TYPE>
    TYPE **grow(TYPE **&array, int n, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE *array)
    {
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 1d array with index from nlo to nhi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create1d_offset(TYPE *&array, int nlo, int nhi, const char *name)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * (nhi-nlo+1);
      array = (TYPE *) smalloc(nbytes,name);
      array -= nlo;
      return array;
    }

  template <typename TYPE>
    TYPE **create1d_offset(TYPE **&array, int nlo, int nhi, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 1d array with index offset
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy1d_offset(TYPE *array, int offset)
    {
      if (array) sfree(&array[offset]);
    }

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1;
      array = (TYPE **) smalloc(nbytes,name);

      bigint n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2;
      }
      return array;
    }

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **grow(TYPE **&array, int n1, int n2, const char *name)
    {
      if (array == NULL) return create(array,n1,n2,name);

      bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
      TYPE *data = (TYPE *) srealloc(array[0],nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1;
      array = (TYPE **) srealloc(array,nbytes,name);

      bigint n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2;
      }
      return array;
    }

  template <typename TYPE>
    TYPE ***grow(TYPE ***&array, int n1, int n2, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 2d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE **array)
    {
      if (array == NULL) return;
      sfree(array[0]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 2d array with a ragged 2nd dimension
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create_ragged(TYPE **&array, int n1, int *n2, const char *name)
    {
      bigint n2sum = 0;
      for (int i = 0; i < n1; i++) n2sum += n2[i];

      bigint nbytes = ((bigint) sizeof(TYPE)) * n2sum;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1;
      array = (TYPE **) smalloc(nbytes,name);

      bigint n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2[i];
      }
      return array;
    }

  template <typename TYPE>
    TYPE ***create_ragged(TYPE ***&array, int n1, int *n2, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   create a 2d array with 2nd index from n2lo to n2hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE **create2d_offset(TYPE **&array, int n1, int n2lo, int n2hi,
                           const char *name)
    {
      int n2 = n2hi - n2lo + 1;
      create(array,n1,n2,name);
      for (int i = 0; i < n1; i++) array[i] -= n2lo;
      return array;
    }

  template <typename TYPE>
    TYPE ***create2d_offset(TYPE ***&array, int n1, int n2lo, int n2hi,
                            const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 2d array with 2nd index offset
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy2d_offset(TYPE **array, int offset)
    {
      if (array == NULL) return;
      sfree(&array[0][offset]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 3d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
      TYPE **plane = (TYPE **) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE **)) * n1;
      array = (TYPE ***) smalloc(nbytes,name);

      int i,j;
      bigint m;
      bigint n = 0;
      for (i = 0; i < n1; i++) {
        m = ((bigint) i) * n2;
        array[i] = &plane[m];
        for (j = 0; j < n2; j++) {
          plane[m+j] = &data[n];
          n += n3;
        }
      }
      return array;
    }

  template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d array
   last 2 dims must stay the same
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, const char *name)
    {
      if (array == NULL) return create(array,n1,n2,n3,name);

      bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
      TYPE *data = (TYPE *) srealloc(array[0][0],nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
      TYPE **plane = (TYPE **) srealloc(array[0],nbytes,name);
      nbytes = ((bigint) sizeof(TYPE **)) * n1;
      array = (TYPE ***) srealloc(array,nbytes,name);

      int i,j;
      bigint m;
      bigint n = 0;
      for (i = 0; i < n1; i++) {
        m = ((bigint) i) * n2;
        array[i] = &plane[m];
        for (j = 0; j < n2; j++) {
          plane[m+j] = &data[n];
          n += n3;
        }
      }
      return array;
    }

  template <typename TYPE>
    TYPE ****grow(TYPE ****&array, int n1, int n2, int n3, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 3d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ***array)
    {
      if (array == NULL) return;
      sfree(array[0][0]);
      sfree(array[0]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 3d array with 1st index from n1lo to n1hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi,
                            int n2, int n3, const char *name)
    {
      int n1 = n1hi - n1lo + 1;
      create(array,n1,n2,n3,name);
      array -= n1lo;
      return array;
    }

  template <typename TYPE>
    TYPE ****create3d_offset(TYPE ****&array, int n1lo, int n1hi,
                             int n2, int n3, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   free a 3d array with 1st index offset
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy3d_offset(TYPE ***array, int offset)
    {
      if (array) destroy(&array[offset]);
    }

/* ----------------------------------------------------------------------
   create a 3d array with
   1st index from n1lo to n1hi inclusive,
   2nd index from n2lo to n2hi inclusive,
   3rd index from n3lo to n3hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi,
                            int n2lo, int n2hi, int n3lo, int n3hi,
                            const char *name)
    {
      int n1 = n1hi - n1lo + 1;
      int n2 = n2hi - n2lo + 1;
      int n3 = n3hi - n3lo + 1;
      create(array,n1,n2,n3,name);

      bigint m = ((bigint) n1) * n2;
      for (bigint i = 0; i < m; i++) array[0][i] -= n3lo;
      for (int i = 0; i < n1; i++) array[i] -= n2lo;
      array -= n1lo;
      return array;
    }

  template <typename TYPE>
    TYPE ****create3d_offset(TYPE ****&array, int n1lo, int n1hi,
                             int n2lo, int n2hi, int n3lo, int n3hi,
                             const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   free a 3d array with all 3 indices offset
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy3d_offset(TYPE ***array,
                          int n1_offset, int n2_offset, int n3_offset)
    {
      if (array == NULL) return;
      sfree(&array[n1_offset][n2_offset][n3_offset]);
      sfree(&array[n1_offset][n2_offset]);
      sfree(&array[n1_offset]);
    }

/* ----------------------------------------------------------------------
   create a 4d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4,
                    const char *name)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3*n4;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1*n2*n3;
      TYPE **cube = (TYPE **) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE **)) * n1*n2;
      TYPE ***plane = (TYPE ***) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE ***)) * n1;
      array = (TYPE ****) smalloc(nbytes,name);

      int i,j,k;
      bigint m1,m2; 
      bigint n = 0;
      for (i = 0; i < n1; i++) {
        m2 = ((bigint) i) * n2;
        array[i] = &plane[m2];
        for (j = 0; j < n2; j++) {
          m1 = ((bigint) i) * n2 + j;
          m2 = ((bigint) i) * n2*n3 + j*n3;
          plane[m1] = &cube[m2];
          for (k = 0; k < n3; k++) {
            m1 = ((bigint) i) * n2*n3 + j*n3 + k;
            cube[m1] = &data[n];
            n += n4;
          }
        }
      }
      return array;
    }

  template <typename TYPE>
    TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4,
                     const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 4d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE ****array)
    {
      if (array == NULL) return;
      sfree(array[0][0][0]);
      sfree(array[0][0]);
      sfree(array[0]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   create a 5d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4,
                     int n5, const char *name)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3*n4*n5;
      TYPE *data = (TYPE *) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE *)) * n1*n2*n3*n4;
      TYPE **level4 = (TYPE **) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE **)) * n1*n2*n3;
      TYPE ***level3 = (TYPE ***) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE ***)) * n1*n2;
      TYPE ****level2 = (TYPE ****) smalloc(nbytes,name);
      nbytes = ((bigint) sizeof(TYPE ****)) * n1;
      array = (TYPE *****) smalloc(nbytes,name);

      int i,j,k,l;
      bigint m1,m2; 
      bigint n = 0;
      for (i = 0; i < n1; i++) {
        m2 = ((bigint) i) * n2;
        array[i] = &level2[m2];
        for (j = 0; j < n2; j++) {
          m1 = ((bigint) i) * n2 + j;
          m2 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3;
          level2[m1] = &level3[m2];
          for (k = 0; k < n3; k++) {
            m1 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3 + k;
            m2 = ((bigint) i) * n2*n3*n4 +
              ((bigint) j) * n3*n4 + ((bigint) k) * n4;
            level3[m1] = &level4[m2];
            for (l = 0; l < n4; l++) {
              m1 = ((bigint) i) * n2*n3*n4 +
                ((bigint) j) * n3*n4 + ((bigint) k) * n4 + l;
              level4[m1] = &data[n];
              n += n5;
            }
          }
        }
      }
      return array;
    }

  template <typename TYPE>
    TYPE ******create(TYPE ******&array, int n1, int n2, int n3, int n4,
                      int n5, const char *name)
    {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 5d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE *****array)
    {
      if (array == NULL) return;
      sfree(array[0][0][0][0]);
      sfree(array[0][0][0]);
      sfree(array[0][0]);
      sfree(array[0]);
      sfree(array);
    }

/* ----------------------------------------------------------------------
   memory usage of arrays, including pointers
------------------------------------------------------------------------- */

  template <typename TYPE>
    bigint usage(TYPE *array, int n)
    {
      UNUSED(array);
      bigint bytes = ((bigint) sizeof(TYPE)) * n;
      return bytes;
    }

  template <typename TYPE>
    bigint usage(TYPE **array, int n1, int n2)
    {
      UNUSED(array);
      bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2;
      bytes += ((bigint) sizeof(TYPE *)) * n1;
      return bytes;
    }

  template <typename TYPE>
    bigint usage(TYPE ***array, int n1, int n2, int n3)
    {
      UNUSED(array);
      bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
      bytes += ((bigint) sizeof(TYPE *)) * n1*n2;
      bytes += ((bigint) sizeof(TYPE **)) * n1;
      return bytes;
    }

  template <typename TYPE>
    bigint usage(TYPE ****array, int n1, int n2, int n3, int n4)
    {
      UNUSED(array);
      bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2*n3*n4;
      bytes += ((bigint) sizeof(TYPE *)) * n1*n2*n3;
      bytes += ((bigint) sizeof(TYPE **)) * n1*n2;
      bytes += ((bigint) sizeof(TYPE ***)) * n1;
      return bytes;
    }
};

}

#endif

/* ERROR/WARNING messages:

E: Failed to allocate %ld bytes for array %s

Your LAMMPS simulation has run out of memory.  You need to run a
smaller simulation or on more processors.

E: Failed to reallocate %ld bytes for array %s

Your LAMMPS simulation has run out of memory.  You need to run a
smaller simulation or on more processors.

E: Cannot create/grow a vector/array of pointers for %s

LAMMPS code is making an illegal call to the templated memory
allocaters, to create a vector or array of pointers.

*/
