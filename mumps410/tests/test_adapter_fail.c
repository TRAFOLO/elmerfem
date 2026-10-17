/* Run with 2 MPI ranks and linked against fake_metis_fails.c. Rank 0 calls the adapter, whose METIS call
   fails; rank 1 waits in a barrier as it would inside MUMPS. The adapter must stop the whole run: the job
   has to end quickly with a non-zero exit code, and neither rank may reach the end of main. */
#include <stdio.h>
#include <mpi.h>
#include <metis.h>

void metis_nodend_(idx_t *n, idx_t *xadj, idx_t *adjncy, idx_t *numflag, idx_t *options4, idx_t *perm, idx_t *iperm);

int main(int argc, char **argv)
{
  int rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    idx_t n = 2, xadj[3] = {1, 2, 3}, adjncy[2] = {2, 1}, numflag = 1, opt4[8] = {0}, perm[2], iperm[2];
    metis_nodend_(&n, xadj, adjncy, &numflag, opt4, perm, iperm);
    printf("test_adapter_fail: ERROR adapter returned after a METIS failure\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("test_adapter_fail: ERROR rank %d reached the end\n", rank);
  fflush(stdout);
  MPI_Finalize();
  return 0;
}
