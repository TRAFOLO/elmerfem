# Vendored license texts

`license_audit.sh` copies these into the bundle's `licenses/`. They live here because neither
the Elmer sources nor the MSYS2 packages provide them.

## Apache-2.0.txt

The Apache License 2.0, verbatim: sha256 `cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30`,
the canonical apache.org text. Copied from the `LICENSE` file of the importlib_metadata 8.4.0
wheel, line endings normalised to LF. Needed for METIS 5.1.0, compiled into `ElmerGrid.exe`:
the fork's `elmergrid/src/metis-5.1.0/LICENSE.txt` only refers to the Apache License, and
Apache-2.0 section 4(a) requires giving recipients a copy of it.

## GPL-3.0.txt, LGPL-3.0.txt

The GNU GPL version 3 and LGPL version 3, verbatim from https://www.gnu.org/licenses/gpl-3.0.txt
(35 149 bytes, sha256 `3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986`) and
https://www.gnu.org/licenses/lgpl-3.0.txt (7 652 bytes, sha256
`e3a994d82e644b03a792a930f574002658412f62407f5fee083f2555c5f23118`), downloaded 2026-10-09.
Shipped when MUMPS is linked: `libelmersolver.dll` then contains Apache-2.0 METIS, which is
compatible with version 3 of the GNU licenses only, so the bundle distributes that library
under LGPL-3.0-or-later and the GPL solver modules under GPL-3.0-or-later (decision 2026-10-09).

## LICENSE.ARPACK.txt

The Rice University BSD license for ARPACK and P_ARPACK. Taken from the arpack-ng `COPYING`
that SciPy 1.16.2 ships (`scipy/sparse/linalg/_eigen/arpack/COPYING`), line endings normalised
to LF, with only the arpack-ng additions removed (the rename note and the Scilab Enterprises,
Octave patch and gentoo patch copyright lines): `libarpack.dll` and `libparpack.dll` are built
from the original Rice code in the fork's `mathlibs/src/`, which carries no license file.

## LICENSE.msmpi.txt

The MIT License of Microsoft MPI (https://github.com/microsoft/Microsoft-MPI) with a header
naming `libmsmpifec.dll`. The MSYS2 msmpi package installs no license file, so the audit
stages this text for `libmsmpifec.dll`.
