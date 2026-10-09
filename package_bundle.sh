#!/bin/bash
# Package + deploy the validated install as the app bundle.
# REWORKED 2026-07-09. Why: the old script (a) dropped only to Trafolo_third_party_dev
# while the live deploy target moved to Trafolo_third_party (production), (b) had no
# retry on Compress-Archive (AV file-lock races produced truncated/no zips twice),
# and (c) had NO prune guard - on 2026-07-03..06 an unpruned 144-module install was
# zipped and deployed to production unnoticed. Limitation: the prune guard checks
# module count / marker files, not full manifest equality.
#
# Usage:  bash package_bundle.sh            # zip + drop DATED folder into _dev
#         PROD=1 bash package_bundle.sh     # additionally replace the PRODUCTION bundle
#                                           # (backup-renamed; CLOSE THE TRAFOLO APP first)
# Env overrides: ELMER_INSTALL, ELMER_SRC, CODE_DIR
set -e
SRC="${ELMER_SRC:-$(cd "$(dirname "$0")" && pwd)}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
CODE="${CODE_DIR:-$SRC/..}"
# Bundle folder/zip name. Must match the app literals (Choke auto_install_dependencies.py +
# thirdpartydownloaddlg.py), the deployed folder, and the ElmerSetup.zip top-level dir on the
# distribution server. Corrected 2026-08-25: the canonical name is
# ElmerFEM-nogui-mpi-Windows-AMD64 — the 2026-07-09 rename to ElmerFEM-TRbuild-Windows-AMD64
# was never adopted by the app and is retired. Verified 2026-08-26 against the actual source:
# Choke hardcodes this exact string at thirdpartydownloaddlg.py:189
# (`os.path.join(third_party_path, "ElmerFEM-nogui-mpi-Windows-AMD64", "bin")`).
NAME="${BUNDLE_NAME:-ElmerFEM-nogui-mpi-Windows-AMD64}"
# Top-level folder INSIDE the zip. Decoupled from NAME (2026-08-25) because CI
# stamps the zip name with the build date, but the app expects a fixed folder
# name when it unpacks the archive. Defaults to NAME so local runs are unchanged.
FOLDER="${FOLDER_NAME:-$NAME}"
# The app's managed third-party dirs under the CURRENT user's %LOCALAPPDATA%
# (was hardcoded to Juris's profile).
LOCAL_APPDATA="$(cygpath -u "${LOCALAPPDATA:-$USERPROFILE\\AppData\\Local}")"
DEV_DROP="$LOCAL_APPDATA/Trafolo_third_party_dev"
PROD_DROP="$LOCAL_APPDATA/Trafolo_third_party"

# ---- prune guard: refuse to package an unpruned install ----
MODCOUNT=$(ls "$INSTALL/share/elmersolver/lib/"*.dll 2>/dev/null | wc -l)
if [ "$MODCOUNT" -ne 18 ] || [ -f "$INSTALL/bin/ViewFactors.exe" ] || [ -d "$INSTALL/lib" ]; then
  echo "ERROR: install is NOT pruned (modules=$MODCOUNT, want 18; ViewFactors/lib check failed)."
  echo "Run prune_install.sh + license_audit.sh + run_gates.sh first."
  exit 1
fi
grep -q "UMFPACK 4.4" "$INSTALL/SOURCE.txt" || { echo "ERROR: SOURCE.txt stale (no UMFPACK) - run license_audit.sh"; exit 1; }
[ -f "$INSTALL/redist/msmpisetup.exe" ] || { echo "ERROR: redist/msmpisetup.exe missing - Choke's thirdpartydownloaddlg.py requires it to install the MS-MPI prerequisite (see README-TRAFOLO.md); rebuild with BUNDLE_MSMPI_REDIST=ON"; exit 1; }

# ---- bundle date from the solver banner (build date, not packaging date) ----
BDATE=$(env -i PATH="$(cygpath -w "$INSTALL/bin");C:\\Windows\\System32" "$INSTALL/bin/ElmerSolver.exe" 2>&1 | grep -oE "Compiled: [0-9-]+" | grep -oE "[0-9-]+$")
[ -n "$BDATE" ] && echo "bundle build date: $BDATE" || { echo "ERROR: could not read solver banner"; exit 1; }

# ---- stage + zip (retry: AV scanners transiently lock freshly copied files) ----
STAGE="$CODE/$FOLDER"
rm -rf "$STAGE"; cp -r "$INSTALL" "$STAGE"
rm -f "$CODE/$NAME.zip"
ok=0
for i in 1 2 3; do
  powershell.exe -NoProfile -Command "Compress-Archive -Path '$(cygpath -w "$STAGE")' -DestinationPath '$(cygpath -w "$CODE/$NAME.zip")' -Force" 2>/dev/null && { ok=1; break; }
  echo "zip attempt $i failed (likely AV lock) - retrying in 5s"; sleep 5
done
[ $ok = 1 ] || { echo "ERROR: zip failed after 3 attempts"; exit 1; }
ls -lh "$CODE/$NAME.zip" | awk '{print "zip:", $5}'

# ---- dev drop (dated folder, side-by-side) ----
mkdir -p "$DEV_DROP"
rm -rf "$DEV_DROP/$NAME-$BDATE"
cp -r "$STAGE" "$DEV_DROP/$NAME-$BDATE"
echo "dev drop: $DEV_DROP/$NAME-$BDATE"

# ---- production deploy (explicit opt-in; backup-rename the old one) ----
if [ "${PROD:-0}" = "1" ]; then
  mkdir -p "$PROD_DROP"
  if [ -d "$PROD_DROP/$NAME" ]; then
    mv "$PROD_DROP/$NAME" "$PROD_DROP/$NAME.bak.$(date +%Y%m%d_%H%M%S)"
    echo "previous production bundle backed up"
  fi
  cp -r "$STAGE" "$PROD_DROP/$NAME"
  echo "PRODUCTION deployed: $PROD_DROP/$NAME"
fi
rm -rf "$STAGE"

# ---- Choke externals staging (hash-sync defusal, G1) ----
# Next to the zip, outside the source checkout.
EXT="$CODE/choke_externals_staging"
rm -rf "$EXT"; mkdir -p "$EXT"
for m in ProcessFields LoadFields htc_udf; do
  cp "$INSTALL/share/elmersolver/lib/$m.dll" "$EXT/"
  cp "$SRC/fem/src/modules/$m.F90" "$EXT/" 2>/dev/null || true
done
echo "Choke externals staged in $EXT (copy into Choke/externals/Elmer/ + verify hash match)"
