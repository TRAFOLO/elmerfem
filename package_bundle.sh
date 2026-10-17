#!/bin/bash
# Package the validated install as the app bundle.
# Refuses an unpruned or unaudited install (prune guard), retries the zip
# (AV scanners transiently lock freshly copied files). Limitation: the prune
# guard checks module count / marker files, not full manifest equality.
#
# Usage:  bash package_bundle.sh            # zip + drop DATED folder into _dev
#         PROD=1 bash package_bundle.sh     # additionally replace the app's installed
#                                           # bundle (backup-renamed; CLOSE THE APP first)
# Env overrides: ELMER_INSTALL, ELMER_SRC, CODE_DIR, BUNDLE_NAME, FOLDER_NAME
set -e
SRC="${ELMER_SRC:-$(cd "$(dirname "$0")" && pwd)}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
CODE="${CODE_DIR:-$SRC/..}"
# Zip file name.
NAME="${BUNDLE_NAME:-ElmerFEM-nogui-mpi-Windows-AMD64}"
# Top-level folder INSIDE the zip. The TRAFOLO app expects exactly
# ElmerFEM-nogui-mpi-Windows-AMD64 when it unpacks the archive, while release zips
# carry the build date in the file name, so the two are set separately.
# Defaults to NAME.
FOLDER="${FOLDER_NAME:-$NAME}"
# The app's third-party dirs under the current user's %LOCALAPPDATA%.
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
[ -f "$INSTALL/licenses/THIRD_PARTY_NOTICES.md" ] || { echo "ERROR: licenses/THIRD_PARTY_NOTICES.md missing (SOURCE.txt points to it) - run license_audit.sh"; exit 1; }
[ -f "$INSTALL/redist/msmpisetup.exe" ] || { echo "ERROR: redist/msmpisetup.exe missing - the app runs it to install the MS-MPI prerequisite (see README-TRAFOLO.md); rebuild with BUNDLE_MSMPI_REDIST=ON"; exit 1; }

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

# ---- stage the TRAFOLO modules for the app's own copies ----
# The app keeps its own copies of these DLLs and must use ones from this same
# build (module DLLs from a different core crash at init). Next to the zip.
EXT="$CODE/choke_externals_staging"
rm -rf "$EXT"; mkdir -p "$EXT"
for m in ProcessFields LoadFields htc_udf; do
  cp "$INSTALL/share/elmersolver/lib/$m.dll" "$EXT/"
  cp "$SRC/fem/src/modules/$m.F90" "$EXT/" 2>/dev/null || true
done
echo "TRAFOLO modules staged in $EXT"
