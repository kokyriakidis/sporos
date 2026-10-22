#!/usr/bin/env bash
# Seed-index lookup benchmark: sporos's SeedIndex (benchmarks/bench_index.mojo) against minimap2's mm_sketch +
# mm_idx_get and strobealign's find_hits (harness/lookup), on 1 M simulated 150 bp reads against hg38 chr21.
#
#     bash benchmarks/competitors/run_lookups.sh
#
# Each harness builds its index, then seeds every read and looks its seeds up on one thread, timing only the
# read loop with the reads in memory. Writes results/lookups.tsv: the best of $REPS runs of each line.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"

R=$D/hg38_chr21.fa
Q=$BENCH/reads/chr21_150bp_1M.fa
LB=$BENCH/mybench/lookup
RAW=$RES/lookups_raw.txt
TSV=$RES/lookups.tsv
MAX_OCC=${MAX_OCC:-20}
{ echo "== run_lookups.sh"; machine_info; } >> "$RES/machine.txt"
: > "$RAW"

for r in $(seq "$REPS"); do
  for kw in "19 12" "21 11"; do
    set -- $kw
    wait_idle
    "$SBIN/bench_index" "$R" "$Q" "$1" "$2" "$MAX_OCC" 2>&1 | sed "s/^/sporos k=$1 w=$2\t/" >> "$RAW"
    wait_idle
    "$LB/mm2_lookup" "$R" "$Q" "$1" "$2" "$MAX_OCC" 2>&1 | sed "s/^/minimap2 k=$1 w=$2\t/" >> "$RAW"
  done
  wait_idle
  "$LB/sa_bench" "$R" "$Q" 2>&1 | sed "s/^/strobealign defaults\t/" >> "$RAW"
done

# Best of the runs: the highest throughput and the lowest build time of each line.
python3 -I - "$RAW" > "$TSV" <<'PY'
import re, sys
best = {}
for line in open(sys.argv[1]):
    tool, _, text = line.rstrip("\n").partition("\t")
    m = re.match(r"(lookup .*?): ([\d.]+) Mbp/s(?: wall)?(?:, ([\d.]+) (?:hits/read|hits per read|per read))?", text)
    if m:
        key = (tool, m.group(1))
        value = (float(m.group(2)), m.group(3) or "")
        if key not in best or value[0] > best[key][0]:
            best[key] = value
        continue
    m = re.match(r"build threads=(\d+): ([\d.]+) s", text)
    if m:
        key = (tool, f"build, {m.group(1)} threads")
        value = (float(m.group(2)), "")
        if key not in best or value[0] < best[key][0]:
            best[key] = value
print("tool\tmeasure\tvalue\thits_per_read")
for (tool, measure), (value, hits) in best.items():
    print(f"{tool}\t{measure}\t{value}\t{hits}")
PY
cat "$TSV"
echo "run_lookups.sh done"
