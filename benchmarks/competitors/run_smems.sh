#!/usr/bin/env bash
# SMEM and f-SMEM benchmark: sporos smem vs minibwa fastmap, bwa fastmap and ropebwt3 mem, on simulated reads
# against hg38 chr21 (setup_xaros.sh makes the reads).
#
#     bash benchmarks/competitors/run_smems.sh
#
# Writes results/smems.tsv (timings, index builds included) and results/smems_checks.tsv (sporos against each
# tool that lists positions). Every tool reads its index from disk; minibwa and bwa fastmap are
# single-threaded, ropebwt3 is run without positions (`-p 0`, its default), which favours it.
# CHECK_ONLY=1 rewrites only results/smems_checks.tsv: every output again, untimed, ropebwt3 skipped.
# SPOROS_ONLY=1 times sporos alone, keeping the other tools' rows, and checks it against the minibwa and bwa
# outputs the last full run kept.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"

THREADS=${THREADS:-"1 10"}
FS=${FS:-"1 2 5"}
L=${L:-19}
SPOROS=$SBIN/sporos
MINIBWA=$BENCH/competitors/minibwa/minibwa
TSV=$RES/smems.tsv
CHK=$RES/smems_checks.tsv
IDX=$BENCH/idx/smem
P=$OUT/smem
mkdir -p "$IDX" "$P"
if [ "${CHECK_ONLY:-0}" = 1 ]; then
  measure() {  # untimed: run the command, its output to the file given (or nowhere)
    while [ "$1" != "--" ]; do shift; done
    shift
    local out=$1; shift
    if [ "$out" = "-" ]; then "$@" > /dev/null; else "$@" > "$out"; fi
  }
  TSV=/dev/null
  rm -f "$CHK"
fi
cached && { drop_rows "$TSV" 3 '^sporos'; rm -f "$CHK"; }
[ -e "$TSV" ] || printf "reads\tL\ttool\tmode\tthreads\treps\twall_s\tuser_s\tsys_s\tmax_rss_mb\tbusy_cpus\tload1\n" > "$TSV"
[ -f "$CHK" ] || printf "reads\tL\tf\tcheck\tn_sporos\tn_other\tonly_sporos\tonly_other\tother_through_N\tother_unlisted\n" > "$CHK"

{ echo "== run_smems.sh"; machine_info; } >> "$RES/machine.txt"

# The reference beside its indexes: sporos looks for <ref>.sfm.
REF=$IDX/chr21.fa
[ -e "$REF" ] || ln -s "$D/hg38_chr21.fa" "$REF"

for T in $THREADS; do
  measure "$TSV" chr21 - "sporos index" build "$T" -- - "$SPOROS" index -t "$T" "$REF"
  cached || measure "$TSV" chr21 - "minibwa index" build "$T" -- - "$MINIBWA" index -t "$T" "$REF" "$IDX/chr21.mb"
  cached || measure "$TSV" chr21 - "ropebwt3 build" build "$T" -- - "$ENVBIN/ropebwt3" build -t "$T" -do "$IDX/chr21.fmd" "$REF"
done
cached || measure "$TSV" chr21 - "bwa index" build 1 -- - "$ENVBIN/bwa" index -p "$IDX/chr21.bwa" "$REF"

for set in 150bp:chr21_150bp_1M 10kb:chr21_10kb_4k; do
  name=${set%%:*}
  Q=$BENCH/reads/${set#*:}.fa
  for f in $FS; do
    mode=$([ "$f" = 1 ] && echo SMEM || echo "f-SMEM f=$f")
    for T in $THREADS; do
      measure "$TSV" "$name" "$L" sporos "$mode" "$T" -- "$P/sporos_${name}_f${f}_T$T.txt" \
        "$SPOROS" smem -l "$L" -s "$f" -t "$T" "$REF" "$Q"
    done
    # minibwa's batched search (-b 64) is its fastest.
    cached || measure "$TSV" "$name" "$L" minibwa "$mode" 1 -- "$P/minibwa_${name}_f$f.txt" \
      "$MINIBWA" fastmap -l "$L" -s "$f" -b 64 "$IDX/chr21.mb" "$Q"
    [ -s "$P/minibwa_${name}_f$f.txt" ] &&
    python3 -I "$HERE/smem_diff.py" "$REF" "$P/sporos_${name}_f${f}_T1.txt" "$P/minibwa_${name}_f$f.txt" |
      awk -v OFS='\t' -v r="$name" -v l="$L" -v f="$f" '{print r, l, f, "sporos vs minibwa", $0}' | tee -a "$CHK"
    if [ "$f" = 1 ]; then
      cached || measure "$TSV" "$name" "$L" "bwa fastmap" "$mode" 1 -- "$P/bwa_${name}.txt" \
        "$ENVBIN/bwa" fastmap -l "$L" "$IDX/chr21.bwa" "$Q"
      [ -s "$P/bwa_${name}.txt" ] &&
      python3 -I "$HERE/smem_diff.py" "$REF" "$P/sporos_${name}_f1_T1.txt" "$P/bwa_${name}.txt" |
        awk -v OFS='\t' -v r="$name" -v l="$L" '{print r, l, 1, "sporos vs bwa fastmap", $0}' | tee -a "$CHK"
    fi
    [ "${CHECK_ONLY:-0}" = 1 ] || cached || for T in $THREADS; do
      measure "$TSV" "$name" "$L" "ropebwt3 mem" "$mode (no positions)" "$T" -- - \
        "$ENVBIN/ropebwt3" mem -t "$T" -l "$L" -c "$f" "$IDX/chr21.fmd" "$Q"
    done
    # Threads change nothing in sporos's output.
    if cmp -s "$P/sporos_${name}_f${f}_T1.txt" "$P/sporos_${name}_f${f}_T10.txt"; then same=PASS; else same=FAIL; fi
    printf "%s\t%s\t%s\tsporos -t10=-t1\t-\t-\t-\t-\t%s\t-\n" "$name" "$L" "$f" "$same" | tee -a "$CHK"
  done
  # minibwa's and bwa's outputs stay for SPOROS_ONLY=1.
  rm -f "$P"/sporos_"$name"_*.txt
done
echo "run_smems.sh done"
