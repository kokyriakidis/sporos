#!/usr/bin/env bash
# Frequency-capped MEM benchmark: sporos mem --max-occ against Mumemto's -f (per-genome cap) and -F (total
# cap) on the genome pairs, both strands, L = 50.
#
#     bash benchmarks/competitors/run_capped.sh
#     PAIRS=ecoli CAPS=3 bash benchmarks/competitors/run_capped.sh
#
# Writes results/capped.tsv (timings) and results/capped_checks.tsv (forward-strand overlap; information
# only). The definitions differ: Mumemto finds multi-MEMs of the two genomes as one collection, counting
# copies on both strands of each, and prints each string once with all its copies; sporos finds the MEMs of
# the query against the reference and counts a string's copies in the reference (and, with -mum, in the query
# sequence). Mumemto -f N pairs with sporos -mum --max-occ N, and -F N with -maxmatch --max-occ N.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"

PAIRS=${PAIRS:-"ecoli chr21"}
CAPS=${CAPS:-"3 10"}
THREADS=${THREADS:-"1 10"}
L=${L:-50}
SPOROS=$SBIN/sporos
TSV=$RES/capped.tsv
CHK=$RES/capped_checks.tsv
MUMEMTO_TIMEOUT=${MUMEMTO_TIMEOUT:-600}
[ -f "$TSV" ] || printf "pair\tL\ttool\tmode\tthreads\treps\twall_s\tuser_s\tsys_s\tmax_rss_mb\tbusy_cpus\tload1\n" > "$TSV"
[ -f "$CHK" ] || printf "pair\tL\tmode\tmumemto_lines\tmumemto_fwd_pairs\tsporos_fwd_pairs\tcommon\n" > "$CHK"

{ echo "== run_capped.sh"; machine_info; } >> "$RES/machine.txt"

for pair in $PAIRS; do
  case $pair in
    ecoli) R=$D/ecoli_k12.fa; Q=$D/ecoli_o157.fa ;;
    chr21) R=$D/hg38_chr21.fa; Q=$D/panTro_chr21_PTRv2.fa ;;
    chr1) R=$D/hg38_chr1.fa; Q=$D/chm13_chr1.fa ;;
  esac
  P=$OUT/capped/$pair
  mkdir -p "$P"
  for N in $CAPS; do
    for cap in "per-genome:-f $N:-mum" "total:-f 0 -F $N:-maxmatch"; do
      kind=${cap%%:*}; rest=${cap#*:}; mm_opts=${rest%%:*}; sp_mode=${rest#*:}
      mode="$kind cap $N"
      for T in $THREADS; do
        measure "$TSV" "$pair" "$L" sporos "$mode" "$T" -- "$P/sporos_${kind}_$N.txt" \
          "$SPOROS" mem "$sp_mode" -b --max-occ "$N" -l "$L" -t "$T" "$R" "$Q"
      done
      rm -f "$P/mumemto_${kind}_$N.mems"
      measure "$TSV" "$pair" "$L" mumemto "$mode" 1 -- - \
        timeout "$MUMEMTO_TIMEOUT" "$ENVBIN/mumemto" -k 2 $mm_opts -l "$L" -o "$P/mumemto_${kind}_$N" "$R" "$Q"
      if [ -s "$P/mumemto_${kind}_$N.mems" ]; then
        python3 -I "$HERE/capped_overlap.py" "$P/mumemto_${kind}_$N.mems" "$P/sporos_${kind}_$N.txt" |
          awk -v OFS='\t' -v p="$pair" -v l="$L" -v m="$mode" '{print p, l, m, $0}' | tee -a "$CHK"
      else
        printf "%s\t%s\t%s\tDNF\tNA\tNA\tNA\n" "$pair" "$L" "$mode" | tee -a "$CHK"
      fi
    done
  done
  rm -f "$P"/*.txt "$P"/*.mems
done
echo "run_capped.sh done"
