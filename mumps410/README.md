# MUMPS 4.10.0 for the TRAFOLO Elmer bundle

`build_mumps410.sh` builds the public-domain MUMPS 4.10.0 (arithmetics s, c, d, z) that
`build_msys2.sh` links into `libelmersolver` when `MUMPS_PREFIX` is set. Run it in the MSYS2
UCRT64 shell with the Microsoft MPI runtime installed; the header lists the paths and packages.

## Why 4.10.0

MUMPS 5.0 and later are distributed under CeCILL-C, whose only declared compatibility is with
CeCILL, not the GPL. MUMPS 4.10.0 (May 2011) is the last release before that change; its
`LICENSE` file declares it public domain. The build installs that file as
`LICENSE.MUMPS-4.10.0.txt` next to the libraries.

## Source

`mumps_4.10.0.dfsg.orig.tar.gz` of Debian (upstream MUMPS 4.10.0 without the user guide),
1 785 786 bytes, SHA-256 `c76339bba516b96a3021af93d9a31b0fbf5a68cfcd02c9578d665ba8018e4b11`.
The script downloads it from the hash-addressed archive
https://snapshot.debian.org/file/5971cd9ccd8c2d789a222140ca6a3b71d6f5229d when it is missing and
always checks the hash.

## What the build changes

- **No PORD.** The orderings are METIS and the MUMPS built-in ones (AMD, AMF, QAMD). An empty
  `libpord.a` only satisfies Elmer's `FindMumps`.
- **METIS 5.1 through `metis4_on_metis5.c`.** MUMPS 4.10.0 calls the METIS 4 entry points; the
  adapter implements them on the METIS 5 API (MSYS2 `mingw-w64-ucrt-x86_64-metis`, Apache-2.0,
  linked statically by `build_msys2.sh`). A METIS failure stops the whole MPI job; there is no
  silent fallback ordering. `tests/` checks both: the result equals a direct `METIS_NodeND` call,
  and a failing METIS ends a 2-rank job with a non-zero exit code.
- **Argument kind of `MUMPS_731`.** `?mumps_part5.F` passes a default `INTEGER` where an
  `INTEGER(8)` is expected; the call now converts it. The build fails if gfortran still reports an
  `INTEGER(4)`/`INTEGER(8)` mismatch.

The adapter, its tests and the source fix are TRAFOLO's, LGPL-2.1-or-later.

## Checks

The script stops at the first failure and prints `PASS` at the end. It checks the archive hash,
that no PORD object is compiled, that every arithmetic calls the METIS ordering, that no ParMETIS
or PORD symbol is present, and records the inputs in `BUILD-INFO.txt` of the install prefix.
