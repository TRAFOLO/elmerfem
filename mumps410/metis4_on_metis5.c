/* METIS 4 ordering entry points used by MUMPS 4.10.0 (Fortran calls METIS_NODEND / METIS_NODEWND),
   implemented on top of the METIS 5 API (METIS 5.1, Apache-2.0), so that the public-domain MUMPS 4.10.0
   can use METIS orderings without the non-commercial METIS 4 library. Written for this purpose; not
   derived from METIS 4 source. Licence of this file: LGPL-2.1-or-later (TRAFOLO).

   MUMPS 4.10.0 calls (src/?mumps_part2.F, ?mumps_part3.F, compiled with -DAdd_):
     METIS_NODEND (N, XADJ, ADJNCY, NUMFLAG, OPTIONS, PERM, IPERM)
     METIS_NODEWND(N, XADJ, ADJNCY, VWGT, NUMFLAG, OPTIONS, PERM, IPERM)
   with default Fortran INTEGER (4 bytes), NUMFLAG = 1 and OPTIONS(1) = 0 (METIS defaults). PERM/IPERM
   have the same meaning in METIS 4 and METIS 5 (row i of the permuted matrix is row PERM(i) of A).

   Failure policy: MUMPS cannot receive an error from these calls, and the ordering runs on one rank
   while the other ranks wait inside MPI. Any METIS failure therefore stops the whole run with
   MPI_Abort (or abort() when MPI is not initialized). There is no silent fallback ordering. */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <metis.h>

_Static_assert(sizeof(idx_t) == 4, "MUMPS 4.10.0 passes 4-byte Fortran INTEGERs; METIS must use 32-bit idx_t");

static void stop_run(const char *what, int status, idx_t n)
{
  int initialized = 0;

  fprintf(stderr, "metis4_on_metis5: %s (status %d, n = %d); stopping the run\n", what, status, (int)n);
  fflush(stderr);
  if (MPI_Initialized(&initialized) == MPI_SUCCESS && initialized)
    MPI_Abort(MPI_COMM_WORLD, 1);
  abort();
}

static void nodend(idx_t *n, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *numflag, idx_t *perm, idx_t *iperm)
{
  idx_t options[METIS_NOPTIONS];
  int status;

  if (*numflag != 0 && *numflag != 1)
    stop_run("invalid NUMFLAG from MUMPS", (int)*numflag, *n);
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = *numflag;
  status = METIS_NodeND(n, xadj, adjncy, vwgt, options, perm, iperm);
  if (status != METIS_OK)
    stop_run("METIS_NodeND failed", status, *n);
}

void metis_nodend_(idx_t *n, idx_t *xadj, idx_t *adjncy, idx_t *numflag, idx_t *options4, idx_t *perm, idx_t *iperm)
{
  (void)options4;
  nodend(n, xadj, adjncy, NULL, numflag, perm, iperm);
}

void metis_nodewnd_(idx_t *n, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *numflag, idx_t *options4,
                    idx_t *perm, idx_t *iperm)
{
  (void)options4;
  nodend(n, xadj, adjncy, vwgt, numflag, perm, iperm);
}
