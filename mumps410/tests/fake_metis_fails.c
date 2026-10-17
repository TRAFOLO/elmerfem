/* Test double: a METIS whose ordering always fails with METIS_ERROR_MEMORY. */
#include <metis.h>

int METIS_SetDefaultOptions(idx_t *options)
{
  for (int i = 0; i < METIS_NOPTIONS; i++)
    options[i] = -1;
  return METIS_OK;
}

int METIS_NodeND(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *options, idx_t *perm, idx_t *iperm)
{
  (void)nvtxs; (void)xadj; (void)adjncy; (void)vwgt; (void)options; (void)perm; (void)iperm;
  return METIS_ERROR_MEMORY;
}
