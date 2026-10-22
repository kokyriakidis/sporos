#!/usr/bin/env bash
# Every benchmark the paper reports, in order, into $RES (default results/), after setup_xaros.sh:
#
#     bash benchmarks/competitors/run_all.sh                  # every tool (about 3 h)
#     SPOROS_ONLY=1 bash benchmarks/competitors/run_all.sh    # sporos alone (about 20 min)
#
# SPOROS_ONLY=1 keeps every other tool's rows in $RES from the last full run and checks sporos against the
# outputs that run kept in ~/sporos-bench/out, so a change to sporos is measured and verified without timing
# the other tools again. Run the full suite again when a tool's version, the machine or the data changes.
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export RES=${RES:-$HERE/results} SPOROS_ONLY=${SPOROS_ONLY:-0}
mkdir -p "$RES"
bash "$HERE/setup_xaros.sh" > "$RES/setup.log" 2>&1 || { echo "setup failed: $RES/setup.log" >&2; exit 1; }
for run in run_seeds run_lookups run_capped run_smems run_mems; do
  echo "== $run $(date -Iseconds)" | tee -a "$RES/run.log"
  bash "$HERE/$run.sh" >> "$RES/run.log" 2>&1
done
echo "ALL DONE $(date -Iseconds)" | tee -a "$RES/run.log"
