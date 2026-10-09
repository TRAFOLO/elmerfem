#!/bin/bash
# ONE-COMMAND driver: source -> validated, packaged TRAFOLO Elmer bundle.
# Runs: build -> deploy -> prune -> audit -> gates -> package, stopping on the
# first failure. Why: skipping a stage is never safe (ninja install silently
# UN-prunes); this driver + package_bundle.sh's prune guard prevent packaging an
# unpruned install. Limitation: Windows/MSYS2-UCRT64 only, by design.
#
# Usage (from any shell - re-execs itself under the UCRT64 login shell):
#   bash build_all.sh                  # full pipeline
#   PROD=1 bash build_all.sh           # + replace the app's installed bundle (close the app!)
#   FROM=gates bash build_all.sh       # resume from: build|deploy|prune|audit|gates|package
#   STOP_AFTER=gates bash build_all.sh # stop early (e.g. validate without packaging)
#   DRYRUN=1 bash build_all.sh         # print the stage plan, run nothing
# Paths (env overrides; defaults are this checkout and its siblings):
#   ELMER_SRC     this repo                      (default: directory of this script)
#   ELMER_BUILD   build dir, WIPED on each build (default: ../elmer-build-win)
#   ELMER_INSTALL install tree, pruned in place  (default: ../elmer-install-win)
#   ELMER_GATES   validation scratch             (default: ../elmer-gates)
#   ELMER_VAL     gate input cases               (default: trafolo_bundle/validation)
#   CODE_DIR      where the zip lands            (default: ..)
# See README-TRAFOLO.md "Building the distributed bundle" for the full procedure.
set -e
# Control flags travel as ARGV through the re-exec: env vars are not reliably
# preserved into the MSYS2 login shell. argv survives everything.
FROM="${1:-${FROM:-build}}"
STOP_AFTER="${2:-${STOP_AFTER:-package}}"
PROD="${3:-${PROD:-0}}"
DRYRUN="${4:-${DRYRUN:-0}}"
if [ "$MSYSTEM" != "UCRT64" ]; then
  SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
  exec /c/msys64/usr/bin/env MSYSTEM=UCRT64 CHERE_INVOKING=1 \
       /c/msys64/usr/bin/bash -l "$SELF" "$FROM" "$STOP_AFTER" "$PROD" "$DRYRUN"
fi
export PROD
HERE="$(cd "$(dirname "$0")" && pwd)"
ORDER="build deploy prune audit gates package"
case " $ORDER " in *" $FROM "*) ;; *) echo "ERROR: FROM=$FROM is not one of: $ORDER"; exit 1;; esac
case " $ORDER " in *" $STOP_AFTER "*) ;; *) echo "ERROR: STOP_AFTER=$STOP_AFTER is not one of: $ORDER"; exit 1;; esac
T0=$(date +%s)
# Every script that turns this source into the distributed zip lives in this repo
# (GPL-2 s.3; decision 2026-10-09). Pin the dirs here so every stage sees the same ones.
export ELMER_SRC="${ELMER_SRC:-$HERE}"
export ELMER_BUILD="${ELMER_BUILD:-$ELMER_SRC/../elmer-build-win}"
export ELMER_INSTALL="${ELMER_INSTALL:-$ELMER_SRC/../elmer-install-win}"
export ELMER_GATES="${ELMER_GATES:-$ELMER_SRC/../elmer-gates}"
export ELMER_VAL="${ELMER_VAL:-$HERE/trafolo_bundle/validation}"
export CODE_DIR="${CODE_DIR:-$ELMER_SRC/..}"

stage_build() {
  # explicit &&: set -e is suspended inside functions called from `if !` below
  bash "$HERE/build_msys2.sh" && ( cd "$ELMER_BUILD" && ninja )
}
stage_deploy()  { bash "$HERE/deploy_msys2.sh"; }
stage_prune()   { bash "$HERE/prune_install.sh"; }
stage_audit()   { bash "$HERE/license_audit.sh"; }
stage_gates()   { bash "$HERE/run_gates.sh"; }
stage_package() { bash "$HERE/package_bundle.sh"; }

started=0
for s in $ORDER; do
  if [ "$s" = "$FROM" ]; then started=1; fi
  if [ "$started" = "0" ]; then
    echo "== [$s] skipped (FROM=$FROM)"
  else
    if [ "$DRYRUN" = "1" ]; then
      echo "== [$s] would run (DRYRUN)"
    else
      t=$(date +%s)
      echo ""
      echo "==================== [$s] started $(date +%H:%M:%S) ===================="
      if ! "stage_$s"; then
        echo "!!!! [$s] FAILED after $(( $(date +%s)-t ))s - fix, then re-run: FROM=$s bash build_all.sh"
        exit 1
      fi
      echo "==================== [$s] OK ($(( $(date +%s)-t ))s) ===================="
    fi
    if [ "$s" = "$STOP_AFTER" ]; then
      echo "STOP_AFTER=$STOP_AFTER reached - done ($(( ($(date +%s)-T0)/60 )) min total)"
      exit 0
    fi
  fi
done

echo ""
echo "ALL STAGES OK in $(( ($(date +%s)-T0)/60 )) min."
