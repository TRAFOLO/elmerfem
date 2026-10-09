/* Checks metis4_on_metis5.c against the real METIS 5.1: on a 2D grid graph, both entry points must return
   valid, mutually inverse permutations in the requested numbering, identical to a direct METIS_NodeND call. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <metis.h>

void metis_nodend_(idx_t *n, idx_t *xadj, idx_t *adjncy, idx_t *numflag, idx_t *options4, idx_t *perm, idx_t *iperm);
void metis_nodewnd_(idx_t *n, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *numflag, idx_t *options4,
                    idx_t *perm, idx_t *iperm);

enum { NX = 30, NY = 20, N = NX * NY };

static idx_t xadj[N + 1], adjncy[4 * N], vwgt[N];

static int build_grid(idx_t base)
{
  idx_t k = 0;
  for (int j = 0; j < NY; j++)
    for (int i = 0; i < NX; i++) {
      int v = j * NX + i;
      xadj[v] = k + base;
      if (i > 0) adjncy[k++] = v - 1 + base;
      if (i < NX - 1) adjncy[k++] = v + 1 + base;
      if (j > 0) adjncy[k++] = v - NX + base;
      if (j < NY - 1) adjncy[k++] = v + NX + base;
      vwgt[v] = 1 + (v % 3);
    }
  xadj[N] = k + base;
  return (int)k;
}

static int check_perm(const char *what, const idx_t *perm, const idx_t *iperm, idx_t base)
{
  static char seen[N];
  memset(seen, 0, sizeof seen);
  for (int i = 0; i < N; i++) {
    idx_t p = perm[i] - base, q = iperm[i] - base;
    if (p < 0 || p >= N || q < 0 || q >= N || seen[p]) {
      printf("FAIL %s: not a permutation at %d\n", what, i);
      return 1;
    }
    seen[p] = 1;
  }
  for (int i = 0; i < N; i++)
    if (perm[iperm[i] - base] - base != i) {
      printf("FAIL %s: perm and iperm are not inverse at %d\n", what, i);
      return 1;
    }
  return 0;
}

static int run(idx_t base)
{
  idx_t n = N, numflag = base, opt4[8] = {0}, options[METIS_NOPTIONS];
  static idx_t p1[N], ip1[N], p2[N], ip2[N], p3[N], ip3[N];
  int bad = 0;

  build_grid(base);
  metis_nodend_(&n, xadj, adjncy, &numflag, opt4, p1, ip1);
  bad |= check_perm("metis_nodend_", p1, ip1, base);
  metis_nodewnd_(&n, xadj, adjncy, vwgt, &numflag, opt4, p2, ip2);
  bad |= check_perm("metis_nodewnd_", p2, ip2, base);
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = base;
  if (METIS_NodeND(&n, xadj, adjncy, NULL, options, p3, ip3) != METIS_OK) {
    printf("FAIL direct METIS_NodeND\n");
    return 1;
  }
  if (memcmp(p1, p3, sizeof p1) || memcmp(ip1, ip3, sizeof ip1)) {
    printf("FAIL numbering %d: adapter result differs from a direct METIS_NodeND call\n", (int)base);
    bad = 1;
  }
  return bad;
}

int main(void)
{
  int bad = run(1) | run(0);
  printf(bad ? "test_adapter_ok: FAILED\n" : "test_adapter_ok: OK\n");
  return bad;
}
