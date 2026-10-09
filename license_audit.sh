#!/bin/bash
# Phase 4 -- license audit gate: every shipped binary's import closure must be
# in the allowlist; forbidden names fail the build. Also stages license texts,
# licenses/THIRD_PARTY_NOTICES.md and SOURCE.txt into the install root.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${ELMER_SRC:-$HERE}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
# Direct solvers: Elmer's VENDORED UMFPACK 4.4 (in-tree under $SRC/umfpack) and, when
# linked, the public-domain MUMPS 4.10.0 of mumps410/build_mumps410.sh. MUMPS 5.x
# (CeCILL-C, GPL-incompatible) never ships: a mumps DLL import fails below, and a
# statically linked MUMPS must be the pinned 4.10.0 (its BUILD-INFO.txt, below).
# git may be absent from the MSYS2 login shell; fall back to Git-for-Windows.
# safe.directory avoids dubious-ownership refusals when the checkout belongs to
# another Windows user.
GIT="$(command -v git || echo '/c/Program Files/Git/cmd/git.exe')"
git() { "$GIT" -c safe.directory='*' "$@"; }
# Prefer the commit the BINARY was built from (banner "Rev:") over the checkout
# HEAD: after doc-only commits, HEAD drifts past the build commit and SOURCE.txt
# would point at a descendant. Fallback to HEAD (still a valid superset pointer).
BANNER=$(env -i PATH="$(cygpath -w "$INSTALL/bin");C:\\Windows\\System32" "$INSTALL/bin/ElmerSolver.exe" 2>/dev/null </dev/null)
BANNER_REV=$(echo "$BANNER" | grep -oE "Rev: [0-9a-f]+" | awk '{print $2}')
# MUMPS is linked statically, so no DLL import shows it: the banner does.
MUMPS_LINKED=0
echo "$BANNER" | grep -q "MUMPS library linked in" && MUMPS_LINKED=1
if [ -n "$BANNER_REV" ] && [ "$BANNER_REV" != "unknown" ]; then
  FORK_COMMIT=$(git -C "$SRC" rev-parse "$BANNER_REV" 2>/dev/null)
fi
[ -z "$FORK_COMMIT" ] && FORK_COMMIT=$(git -C "$SRC" rev-parse HEAD 2>/dev/null)
if [ -z "$FORK_COMMIT" ]; then
  echo "ERROR: could not resolve fork commit (git failed in this shell) -- SOURCE.txt would be incomplete"
  exit 1
fi
echo "SOURCE commit: $FORK_COMMIT (banner Rev: ${BANNER_REV:-n/a})"
FAIL=0

# The staging below deletes licenses/ and SOURCE.txt: never run it in the wrong dir.
cd "$INSTALL" || { echo "ERROR: install dir not found: $INSTALL"; exit 1; }
ALL_IMPORTS=$(for f in bin/*.exe bin/*.dll share/elmersolver/lib/*.dll; do
  objdump -x "$f" 2>/dev/null | grep "DLL Name" | awk '{print $3}'; done | sort -u)

echo "=== aggregated import closure ==="; echo "$ALL_IMPORTS"

# forbidden anywhere in the closure (umfpack-as-DLL forbidden too: only the
# vendored in-tree copy may be used, which links statically and never imports)
for bad in mumps parmetis scotch umfpack suitesparse cholmod fftw readline gsl mkl pardiso hypre zoltan; do
  if echo "$ALL_IMPORTS" | grep -qi "$bad"; then echo "FORBIDDEN IMPORT: $bad"; FAIL=1; fi
done
# msmpi may be imported (system prerequisite) but must not be SHIPPED
if ls bin/msmpi.dll 2>/dev/null; then echo "FORBIDDEN FILE: bin/msmpi.dll must not be zipped"; FAIL=1; fi

# every shipped DLL file must be in the allowlist (catches accidental copies)
# libmsmpifec: since MSYS2 msmpi 10.1.1-23 the MS-MPI Fortran interface is in
# this DLL. MIT (Microsoft-MPI), so it ships with its license text;
# msmpi.dll itself stays a system prerequisite (checked above).
ALLOW='^(libelmersolver|libmatc|libfhuti|libmpi_stubs|libgcc_s_seh-1|libgfortran-5|libquadmath-0|libwinpthread-1|libgomp-1|libstdc\+\+-6|libopenblas|libscalapack|libarpack|libparpack|zlib1|libmsmpifec)\.dll$'
for f in bin/*.dll; do
  b=$(basename "$f")
  echo "$b" | grep -qE "$ALLOW" || { echo "UNLISTED SHIPPED DLL: $b (classify before shipping)"; FAIL=1; }
done

# --- stage license texts + SOURCE note ---
# pacman names the package behind each runtime DLL below; check before deleting anything.
command -v pacman >/dev/null || { echo "ERROR: pacman not found - run the audit in the MSYS2 UCRT64 shell (build_all.sh does)"; exit 1; }
# From scratch: texts of components no longer shipped, or a stale SOURCE.txt,
# must not survive in a reused install dir.
echo "=== staging licenses/ + SOURCE.txt ==="
rm -rf licenses SOURCE.txt
mkdir -p licenses
cp -r "$SRC/license_texts"/* licenses/ 2>/dev/null
# A missing license text fails the audit; it used to be skipped silently.
copy_required() { # src dst
  cp "$1" "$2" 2>/dev/null || { echo "MISSING LICENSE TEXT: $1 (for $2)"; FAIL=1; }
}
# Third-party code compiled into the Elmer binaries: notices from the fork's tree.
copy_required "$SRC/umfpack/License" licenses/LICENSE.UMFPACK.txt
copy_required "$SRC/fem/LICENSES" licenses/LICENSES.fem.txt
copy_required "$SRC/elmergrid/LICENSES" licenses/LICENSES.elmergrid.txt
copy_required "$SRC/elmergrid/src/metis-5.1.0/LICENSE.txt" licenses/LICENSE.METIS.txt
copy_required "$SRC/contrib/lua-5.1.5/COPYRIGHT" licenses/LICENSE.Lua.txt
# Texts the fork does not carry (origin: trafolo_bundle/licenses/README.md).
copy_required "$HERE/trafolo_bundle/licenses/Apache-2.0.txt" licenses/Apache-2.0.txt
copy_required "$HERE/trafolo_bundle/licenses/LICENSE.ARPACK.txt" licenses/LICENSE.ARPACK.txt

# --- runtime DLLs from MSYS2: package, version and license texts of each ---
# deploy_msys2.sh copies with cp -n, so a reused install dir can hold stale DLLs:
# only a byte-identical file shows which package version the notice must name.
PFX="${MSYSTEM_PREFIX:-/ucrt64}"
ELMER_BUILT='^(libelmersolver|libmatc|libfhuti|libmpi_stubs|libarpack|libparpack)\.dll$'
# cmp is in diffutils, which the MSYS2 base install lacks; sha256sum is in coreutils.
same_bytes() { [ "$(sha256sum < "$1")" = "$(sha256sum < "$2")" ]; }
stage_pkg_texts() { # pkg short -> stages into licenses/third_party/<short>/, prints the names
  pacman -Qlq "$1" | grep '/share/licenses/' | while IFS= read -r l; do
    [ -f "$l" ] && mkdir -p "licenses/third_party/$2" && cp "$l" "licenses/third_party/$2/" \
      && basename "$l"
  done
}
declare -A PKG_TEXTS
RT_ROWS=""
# Exact source archives of the MSYS2 packages shipped or linked: the GPL/LGPL corresponding
# source of the runtime DLLs and the static METIS. release.yml mirrors them next to the bundle.
MSYS2_SOURCES="https://repo.msys2.org/mingw/sources"
declare -A SRC_ARCHIVES
add_source_archive() { # pkg version
  local base
  base=$(sed -n '/^%BASE%/{n;p}' "/var/lib/pacman/local/$1-$2/desc" 2>/dev/null)
  [ -n "$base" ] || { echo "ERROR: no source package name (%BASE%) for $1 $2 in the pacman database"; FAIL=1; return; }
  SRC_ARCHIVES["$base-$2.src.tar.zst"]="$MSYS2_SOURCES/$base-$2.src.tar.zst"
}
echo "=== MSYS2 runtime DLLs in bin/ ==="
for f in bin/*.dll; do
  b=$(basename "$f")
  echo "$b" | grep -qE "$ELMER_BUILT" && continue
  if [ ! -f "$PFX/bin/$b" ] || ! same_bytes "$PFX/bin/$b" "$f"; then
    [ -f "$PFX/bin/$b" ] && why="differs from $PFX/bin/$b" || why="has no $PFX/bin/$b"
    echo "FOREIGN/STALE DLL: bin/$b $why - the shipped DLL does not come from the installed MSYS2 package, so its notice would name the wrong version (redeploy into a fresh install dir)"
    FAIL=1; continue
  fi
  owner=$(LC_ALL=C pacman -Qo "$PFX/bin/$b" 2>/dev/null | sed -n 's/^.* is owned by //p')
  [ -n "$owner" ] || { echo "ERROR: no MSYS2 package owns $PFX/bin/$b"; FAIL=1; continue; }
  pkg=${owner% *}; ver=${owner##* }; short=${pkg#mingw-w64-ucrt-x86_64-}
  if [ -z "${PKG_TEXTS[$pkg]+seen}" ]; then
    PKG_TEXTS[$pkg]=$(stage_pkg_texts "$pkg" "$short")
    # No license file in the package (today: msmpi): use the vendored text, which
    # also stays at licenses/LICENSE.<short>.txt as before.
    V="$HERE/trafolo_bundle/licenses/LICENSE.$short.txt"
    if [ -z "${PKG_TEXTS[$pkg]}" ] && [ -f "$V" ] && mkdir -p "licenses/third_party/$short" \
       && cp "$V" "licenses/third_party/$short/" && cp "$V" licenses/; then
      PKG_TEXTS[$pkg]="LICENSE.$short.txt"
    fi
  fi
  [ -n "${PKG_TEXTS[$pkg]}" ] || { echo "MISSING LICENSE TEXT: no license text for $b ($pkg ships none; vendor it as trafolo_bundle/licenses/LICENSE.$short.txt)"; FAIL=1; }
  spdx=$(LC_ALL=C pacman -Qi "$pkg" | sed -n 's/^Licenses *: *//p' | sed 's/spdx://g; s/   */, /g')
  texts=""; for t in ${PKG_TEXTS[$pkg]}; do texts+="${texts:+, }\`$t\`"; done
  RT_ROWS+="| \`bin/$b\` | $pkg $ver | $spdx | \`third_party/$short/\`: $texts | https://packages.msys2.org/packages/$pkg |"$'\n'
  add_source_archive "$pkg" "$ver"
  echo "$b: $pkg $ver ($spdx)"
done

# --- MUMPS 4.10.0 and its METIS ordering, linked statically into libelmersolver ---
# Only the public-domain MUMPS 4.10.0 built by the fork's mumps410/build_mumps410.sh may
# ship; its BUILD-INFO.txt in MUMPS_PREFIX names the pinned source archive and the METIS
# package linked statically for the orderings (Apache-2.0).
MUMPS_SHA256=c76339bba516b96a3021af93d9a31b0fbf5a68cfcd02c9578d665ba8018e4b11
METIS_PKG=mingw-w64-ucrt-x86_64-metis
MUMPS_ROWS=""
if [ "$MUMPS_LINKED" = 1 ]; then
  INFO="${MUMPS_PREFIX:-/nonexistent}/BUILD-INFO.txt"
  if [ ! -s "$INFO" ]; then
    echo "MUMPS: the solver links MUMPS, but $INFO is missing (only MUMPS 4.10.0 from the fork's mumps410/build_mumps410.sh may ship)"; FAIL=1
  elif ! grep -q "^tarball_sha256=$MUMPS_SHA256\$" "$INFO"; then
    echo "MUMPS: $INFO does not name the pinned MUMPS 4.10.0 source archive"; FAIL=1
  else
    copy_required "$MUMPS_PREFIX/LICENSE.MUMPS-4.10.0.txt" licenses/LICENSE.MUMPS-4.10.0.txt
    metis_ver=$(sed -n "s/^metis=$METIS_PKG \([^,]*\),.*/\1/p" "$INFO")
    installed=$(pacman -Q "$METIS_PKG" 2>/dev/null | awk '{print $2}')
    if [ -z "$metis_ver" ] || [ "$metis_ver" != "$installed" ]; then
      echo "MUMPS: METIS linked ($metis_ver) is not the installed $METIS_PKG ($installed): its notice would name the wrong version"; FAIL=1
    fi
    # The MSYS2 package installs no licence file; METIS 5.1.0's own LICENSE.txt (staged above
    # from the fork's ElmerGrid copy of the same release) names the Apache License 2.0.
    metis_texts=$(stage_pkg_texts "$METIS_PKG" metis)
    texts="\`LICENSE.METIS.txt\`"; for t in $metis_texts; do texts+=", \`third_party/metis/$t\`"; done
    MUMPS_ROWS+="| MUMPS 4.10.0 (arithmetics s, c, d, z) | public domain | \`bin/libelmersolver.dll\` | \`LICENSE.MUMPS-4.10.0.txt\` |"$'\n'
    MUMPS_ROWS+="| METIS 5.1.0, MSYS2 $METIS_PKG $metis_ver (MUMPS orderings) | Apache-2.0 | \`bin/libelmersolver.dll\` | $texts, \`Apache-2.0.txt\` |"$'\n'
    MUMPS_ROWS+="| TRAFOLO METIS-4 adapter and MUMPS source fix (fork \`mumps410/\`) | LGPL-2.1-or-later | \`bin/libelmersolver.dll\` | \`LGPL-2.1.txt\` |"$'\n'
    # Apache-2.0 METIS inside libelmersolver is compatible with version 3 of the GNU licenses
    # only: the bundle uses the "or any later version" terms (decision 2026-10-09).
    copy_required "$HERE/trafolo_bundle/licenses/GPL-3.0.txt" licenses/GPL-3.0.txt
    copy_required "$HERE/trafolo_bundle/licenses/LGPL-3.0.txt" licenses/LGPL-3.0.txt
    add_source_archive "$METIS_PKG" "$metis_ver"
    SRC_ARCHIVES["mumps_4.10.0.dfsg.orig.tar.gz"]="https://snapshot.debian.org/file/5971cd9ccd8c2d789a222140ca6a3b71d6f5229d"
    echo "MUMPS 4.10.0 linked (archive SHA-256 pinned), METIS $metis_ver"
  fi
fi

# --- Microsoft MPI redistributable installer ---
# msmpisetup.exe comes under Microsoft's redistributable licence terms, not under the MIT
# License of the Microsoft-MPI source repository. The MS-MPI runtime installed on the build
# machine comes from the same installer (release.yml), so its License folder holds those terms.
MSMPI_LICENSES="${MSMPI_LICENSE_DIR:-/c/Program Files/Microsoft MPI/License}"
if [ -f redist/msmpisetup.exe ]; then
  mkdir -p licenses/third_party/msmpi-redist
  copy_required "$MSMPI_LICENSES/MicrosoftMPI-Redistributable-EULA.rtf" licenses/third_party/msmpi-redist/
  copy_required "$MSMPI_LICENSES/MPI-Redistributables-TPN.txt" licenses/third_party/msmpi-redist/
fi

# --- licenses/THIRD_PARTY_NOTICES.md (SOURCE.txt points here) ---
# The compiled-in table is static and holds only for build_msys2.sh's flags (in-tree
# UMFPACK/AMD, ARPACK/PARPACK and METIS, WITH_LUA=ON): change them together.
NOTICES=licenses/THIRD_PARTY_NOTICES.md
cat > "$NOTICES" <<'EOF'
# Third-party notices

The ElmerFEM solver bundle for TRAFOLO is a Windows build (MSYS2 UCRT64) of the TRAFOLO fork of
Elmer FEM. The Elmer libraries libelmersolver, matc and fhutiter (`bin/libelmersolver.dll`,
`bin/libmatc.dll`, `bin/libfhuti.dll`) are licensed under LGPL-2.1-or-later (`LGPL-2.1.txt`);
ElmerGrid and the GPL solver modules are licensed under GPL-2.0-or-later (`GPL-2.txt`). See
`ElmerLicensePolicy.md`, `LICENSES` and `LICENSES_GPL.txt`. The source code is at
https://github.com/TRAFOLO/elmerfem, at the commit recorded in `SOURCE.txt` in the bundle root.

License file names below are relative to this directory (`licenses/`); `bin/` and `redist/` are
relative to the bundle root.

## Third-party code compiled into the Elmer binaries

| Component | License | Compiled into | License text |
|---|---|---|---|
| UMFPACK 4.4 | LGPL-2.1-or-later | `bin/libelmersolver.dll` | `LICENSE.UMFPACK.txt` |
| AMD 1.1 | AMD License | `bin/libelmersolver.dll` | `LICENSES.fem.txt` (section lic_amd) |
| Lua 5.1.5 | MIT | `bin/libelmersolver.dll` | `LICENSE.Lua.txt` |
| ARPACK, PARPACK | BSD (Rice University) | `bin/libarpack.dll`, `bin/libparpack.dll` | `LICENSE.ARPACK.txt` |
| METIS 5.1.0 | Apache-2.0 | `bin/ElmerGrid.exe` | `LICENSE.METIS.txt`, `Apache-2.0.txt` |

The exception by CSC - IT Center for Science Ltd. that allows linking the GPL-licensed ElmerGrid
with METIS is in `LICENSES.elmergrid.txt`.

The UMFPACK and AMD licenses require user documentation to cite the copyright, the license, the
availability note and "Used by permission.":

```text
UMFPACK Version 4.4, Copyright 1995-2005 by Timothy A. Davis.
All Rights Reserved.
Used by permission.
Availability: http://www.cise.ufl.edu/research/sparse/umfpack

AMD Version 1.1 (Jan. 21, 2004),  Copyright (c) 2004 by Timothy A.
Davis, Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
Used by permission.
Availability: http://www.cise.ufl.edu/research/sparse/amd
```
EOF
if [ -n "$MUMPS_ROWS" ]; then
  {
    printf '\n## MUMPS direct solver, compiled into `bin/libelmersolver.dll`\n\n'
    printf '%s\n' "MUMPS 4.10.0 is the last MUMPS release before MUMPS 5.0 moved to the CeCILL-C license; it is" \
      "public domain. It was built by \`mumps410/build_mumps410.sh\` of the fork from Debian's" \
      "\`mumps_4.10.0.dfsg.orig.tar.gz\` (SHA-256 \`$MUMPS_SHA256\`," \
      "https://snapshot.debian.org/file/5971cd9ccd8c2d789a222140ca6a3b71d6f5229d), with METIS for the" \
      "orderings and without PORD." ""
    printf '| Component | License | Compiled into | License text |\n|---|---|---|---|\n'
    printf '%s' "$MUMPS_ROWS"
    printf '\n## License version of this bundle\n\n'
    printf '%s\n' "\`bin/libelmersolver.dll\` contains METIS (Apache-2.0), which is compatible with version 3 of" \
      "the GNU licenses only. TRAFOLO therefore distributes \`bin/libelmersolver.dll\` under" \
      "LGPL-3.0-or-later (\`LGPL-3.0.txt\`, which adds permissions to \`GPL-3.0.txt\`) and the" \
      "GPL-2.0-or-later Elmer solver modules that run with it under GPL-3.0-or-later (\`GPL-3.0.txt\`)," \
      "as the \"or (at your option) any later version\" terms of their licenses permit. The source files" \
      "keep their original license headers."
  } >> "$NOTICES"
fi
{
  printf '# Exact source archives of the third-party code shipped in or linked into this bundle\n'
  printf '# (MSYS2 source packages: PKGBUILD, patches, upstream tarball). <file> <origin>\n'
  for f in "${!SRC_ARCHIVES[@]}"; do printf '%s %s\n' "$f" "${SRC_ARCHIVES[$f]}"; done | sort
} > licenses/source-packages.txt
if [ -n "${SOURCE_MIRROR_URL:-}" ]; then
  SOURCES_SENTENCE="The exact source archives of these packages and of the third-party code compiled in from MSYS2 or Debian archives (PKGBUILD, patches and upstream tarball), listed in \`source-packages.txt\`, are published with this bundle at $SOURCE_MIRROR_URL and by their origins."
else
  SOURCES_SENTENCE="The exact source archives of these packages and of the third-party code compiled in from MSYS2 or Debian archives (PKGBUILD, patches and upstream tarball) are listed in \`source-packages.txt\` and published by their origins."
fi
{
  printf '\n## Runtime libraries from MSYS2 (UCRT64)\n\n'
  printf '%s\n' "These DLLs in \`bin/\` are unmodified files of the MSYS2 UCRT64 packages listed. $SOURCES_SENTENCE" ""
  printf '| File | Package and version | License (SPDX) | License texts | Source |\n|---|---|---|---|---|\n'
  printf '%s' "$RT_ROWS"
  printf '\n## Microsoft MPI\n\n'
  [ -f bin/libmsmpifec.dll ] && printf '%s\n' \
    "- \`bin/libmsmpifec.dll\`: the MS-MPI Fortran interface, built by MSYS2 from Microsoft MPI." \
    "  MIT License, text in \`LICENSE.msmpi.txt\`."
  [ -f redist/msmpisetup.exe ] && printf '%s\n' \
    "- \`redist/msmpisetup.exe\`: the unmodified Microsoft MPI redistributable installer as published" \
    "  by Microsoft, under Microsoft's redistributable licence terms" \
    "  (\`third_party/msmpi-redist/MicrosoftMPI-Redistributable-EULA.rtf\`, third-party notices in" \
    "  \`third_party/msmpi-redist/MPI-Redistributables-TPN.txt\`). The MIT License of the Microsoft-MPI" \
    "  source repository does not cover this installer."
  printf '%s\n' "- \`msmpi.dll\` itself is not shipped."
} >> "$NOTICES"

if [ -n "$MUMPS_ROWS" ]; then
  MUMPS_NOTE="MUMPS 4.10.0 (direct solver, public domain): compiled into
  bin/libelmersolver.dll by mumps410/build_mumps410.sh of the fork from
  Debian's mumps_4.10.0.dfsg.orig.tar.gz (SHA-256 $MUMPS_SHA256,
  https://snapshot.debian.org/file/5971cd9ccd8c2d789a222140ca6a3b71d6f5229d),
  with METIS 5.1.0 (MSYS2 $METIS_PKG $metis_ver, Apache-2.0) for the
  orderings. MUMPS 5.x (CeCILL-C) is not used. Because Apache-2.0 METIS
  is linked into bin/libelmersolver.dll, that library is distributed
  under LGPL-3.0-or-later and the GPL solver modules under
  GPL-3.0-or-later (licenses/LGPL-3.0.txt, licenses/GPL-3.0.txt)."
else
  MUMPS_NOTE="No MUMPS in this bundle."
fi
cat > SOURCE.txt <<EOF
ElmerFEM solver bundle for TRAFOLO -- source and license notes
==============================================================
Built: $(date -u +%Y-%m-%d)

Elmer sources (GPL-2.0+/LGPL-2.1, see licenses/):
  TRAFOLO fork of ElmerFEM, branch 'trafolo'
  commit $FORK_COMMIT
  Repository: https://github.com/TRAFOLO/elmerfem
  Includes TRAFOLO modules ProcessFields, LoadFields, htc_udf
  (fem/src/modules/, GPL-2.0+). Corresponding sources are available
  from the repository at the commit above, including build_all.sh and
  the scripts it runs to build this bundle.

UMFPACK 4.4 (direct solver; vendored in-tree at umfpack/ of the fork,
  LGPL-2.1-or-later with a user-documentation citation duty):
  part of the corresponding sources above - see licenses/ and
  licenses/THIRD_PARTY_NOTICES.md.

$MUMPS_NOTE

Dynamic third-party runtimes in bin/:
  OpenBLAS (BSD-3), GCC runtimes libgcc/libgfortran/libquadmath
  (GPL with runtime exception), libwinpthread (MIT/BSD), the MS-MPI
  Fortran interface libmsmpifec (MIT, built by MSYS2 from Microsoft-MPI)
  when present, and any further DLLs enumerated by the audit above (see
  licenses/ for their texts).

Third-party source archives: licenses/source-packages.txt lists the
  exact source archive of every MSYS2 package shipped or linked (and of
  MUMPS when linked).${SOURCE_MIRROR_URL:+ They are published with this bundle at
  $SOURCE_MIRROR_URL}

Microsoft MPI runtime (msmpi.dll): not shipped in bin/. The unmodified
  Microsoft installer redist/msmpisetup.exe is included to install it,
  under Microsoft's redistributable licence terms
  (licenses/third_party/msmpi-redist/).
EOF
echo "SOURCE.txt + licenses/ staged."

# Every text SOURCE.txt and THIRD_PARTY_NOTICES.md point to must be there.
for req in THIRD_PARTY_NOTICES.md Apache-2.0.txt LICENSE.METIS.txt LICENSE.Lua.txt \
           LICENSE.ARPACK.txt LICENSE.UMFPACK.txt GPL-2.txt LGPL-2.1.txt; do
  [ -s "licenses/$req" ] || { echo "MISSING FROM licenses/: $req"; FAIL=1; }
done
if [ -f redist/msmpisetup.exe ]; then
  [ -s licenses/third_party/msmpi-redist/MicrosoftMPI-Redistributable-EULA.rtf ] \
    || { echo "MISSING FROM licenses/: third_party/msmpi-redist/MicrosoftMPI-Redistributable-EULA.rtf"; FAIL=1; }
fi
if [ "$MUMPS_LINKED" = 1 ]; then
  for req in LICENSE.MUMPS-4.10.0.txt GPL-3.0.txt LGPL-3.0.txt; do
    [ -s "licenses/$req" ] || { echo "MISSING FROM licenses/: $req"; FAIL=1; }
  done
fi

[ $FAIL = 0 ] && echo "LICENSE AUDIT: CLEAN" || echo "LICENSE AUDIT: FAILED"
exit $FAIL
