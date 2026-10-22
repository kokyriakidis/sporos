#!/usr/bin/env bash
# MEM/MUM benchmark: sporos mem vs MUMmer4, copMEM2, mumemto, bwa fastmap, ropebwt3 mem.
#
#     bash benchmarks/competitors/run_mems.sh            # all pairs
#     PAIRS="ecoli chr21" bash benchmarks/competitors/run_mems.sh
#
# Writes results/mems.tsv (timings) and results/mems_checks.tsv (output equality). SPOROS_ONLY=1 times sporos
# alone, keeping the other tools' rows, and checks it against their sorted outputs from the last full run.
# Best wall time of $REPS runs (a first run over $REPEAT_UNDER s is not repeated).
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"

PAIRS=${PAIRS:-"ecoli chr21 chr1"}
LS=${LS:-"50 100"}
THREADS=${THREADS:-"1 10"}
MODES="maxmatch mum mumreference"
SPOROS=$SBIN/sporos
COPMEM2=$BENCH/competitors/copmem2/copmem2
TSV=$RES/mems.tsv
CHK=$RES/mems_checks.tsv
IDX=$BENCH/idx
# mumemto needs > 300 s on the chr1 pair (> 5x the slowest other MEM tool); cap it.
MUMEMTO_TIMEOUT=${MUMEMTO_TIMEOUT:-300}
mkdir -p "$IDX"
cached && { drop_rows "$TSV" 3 '^sporos$'; rm -f "$CHK"; }
[ -f "$TSV" ] || printf "pair\tL\ttool\tmode\tthreads\treps\twall_s\tuser_s\tsys_s\tmax_rss_mb\tbusy_cpus\tload1\n" > "$TSV"
[ -f "$CHK" ] || printf "pair\tL\tcheck\tn_a\tn_b\tonly_a\tonly_b\tresult\n" > "$CHK"

pair_files() {
  case $1 in
    ecoli) R=$D/ecoli_k12.fa; Q=$D/ecoli_o157.fa ;;
    chr21) R=$D/hg38_chr21.fa; Q=$D/panTro_chr21_PTRv2.fa ;;
    chr1) R=$D/hg38_chr1.fa; Q=$D/chm13_chr1.fa ;;
    *) echo "unknown pair $1" >&2; exit 1 ;;
  esac
}

# MUMmer-format output (3 or 4 columns, '> name [Reverse]' headers) -> sorted set of
# query, strand, ref pos, query pos, length (all references here are single records).
norm() {
  # (mawk 1.3.4 mishandles `$NF == "Reverse"` here, hence index().)
  awk '/^>/ { q = $2; s = index($0, " Reverse") ? "-" : "+"; next }
       NF >= 3 { print q "\t" s "\t" $(NF-2) "\t" $(NF-1) "\t" $NF }' "$1" |
    LC_ALL=C sort -u -S 8G --parallel=8 > "$2"
}

# mumemto .mums (len, 0-based starts, strands) -> forward-strand ref pos, query pos, length.
norm_mumemto() {
  awk -F'\t' '$3=="+,+"{split($2,p,","); print p[1]+1"\t"p[2]+1"\t"$1}' "$1" | LC_ALL=C sort -u > "$2"
}

# check <pair> <L> <name> <set a> <set b> [exact|info]
check() {
  local na nb oa ob res
  na=$(wc -l < "$4"); nb=$(wc -l < "$5")
  oa=$(LC_ALL=C comm -23 "$4" "$5" | wc -l); ob=$(LC_ALL=C comm -13 "$4" "$5" | wc -l)
  if [ "${6:-exact}" = info ]; then res=info; elif [ "$oa" -eq 0 ] && [ "$ob" -eq 0 ]; then res=PASS; else res=FAIL; fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$na" "$nb" "$oa" "$ob" "$res" | tee -a "$CHK"
}

# Replace the timing fields of the last row of a TSV with a DNF note.
mark_dnf() {
  awk -F'\t' -v OFS='\t' -v note="$2" 'NR > 1 { print prev } { prev = $0 }
    END { n = split(prev, f, "\t"); f[6] = 1; f[7] = note; f[8] = f[9] = f[10] = "NA"
          line = f[1]; for (i = 2; i <= n; i++) line = line OFS f[i]; print line }' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}

{ echo "== run_mems.sh"; machine_info; } >> "$RES/machine.txt"

for pair in $PAIRS; do
  pair_files "$pair"
  P=$OUT/$pair
  mkdir -p "$P"

  # FM-index builds, timed once (single thread for bwa; ropebwt3 at each thread count).
  cached || measure "$TSV" "$pair" - "bwa index" build 1 -- - "$ENVBIN/bwa" index -p "$IDX/$pair.bwa" "$R"
  for T in $THREADS; do
    cached || measure "$TSV" "$pair" - "ropebwt3 build" build "$T" -- - "$ENVBIN/ropebwt3" build -t "$T" -do "$IDX/$pair.fmd" "$R"
  done

  for L in $LS; do
    for mode in $MODES; do
      for T in $THREADS; do
        measure "$TSV" "$pair" "$L" sporos "$mode" "$T" -- "$P/sporos_${mode}_L${L}_T$T.txt" \
          "$SPOROS" mem -$mode -l "$L" -n -t "$T" "$R" "$Q"
        norm "$P/sporos_${mode}_L${L}_T$T.txt" "$P/sporos_${mode}_L${L}_T$T.set"
      done
      # MUMmer 4.0.1: -qthreads does not split a single query and -k>1 -threads aborts,
      # so it runs single-threaded.
      cached || { measure "$TSV" "$pair" "$L" mummer "$mode" 1 -- "$P/mummer_${mode}_L$L.txt" \
        "$ENVBIN/mummer" -$mode -n -l "$L" "$R" "$Q"
        norm "$P/mummer_${mode}_L$L.txt" "$P/mummer_${mode}_L$L.set"; }
      check "$pair" "$L" "sporos=mummer -$mode" "$P/sporos_${mode}_L${L}_T1.set" "$P/mummer_${mode}_L$L.set"
      check "$pair" "$L" "sporos -t10=-t1 -$mode" "$P/sporos_${mode}_L${L}_T10.set" "$P/sporos_${mode}_L${L}_T1.set"
    done

    for T in $THREADS; do
      for mf in "" "-mf"; do
        cached || { measure "$TSV" "$pair" "$L" "copmem2${mf:+ $mf}" maxmatch "$T" -- - \
          "$COPMEM2" -q $mf -t "$T" -l "$L" -o "$P/copmem2${mf}_L${L}_T$T.txt" "$R" "$Q"
          norm "$P/copmem2${mf}_L${L}_T$T.txt" "$P/copmem2${mf}_L${L}_T$T.set"; }
        check "$pair" "$L" "sporos=copmem2${mf:+ $mf} -t$T" "$P/sporos_maxmatch_L${L}_T1.set" "$P/copmem2${mf}_L${L}_T$T.set"
      done
    done

    # mumemto: multi-MUMs over both strands (unique among both strands of both inputs);
    # single-threaded; builds its index (PFP + SA/LCP) every run.
    # GNU timeout signals the whole process group (python wrapper and mumemto_exec).
    cached || { rm -f "$P/mumemto_L$L.mums"
      measure "$TSV" "$pair" "$L" mumemto "mum (both strands)" 1 -- - \
        timeout "$MUMEMTO_TIMEOUT" "$ENVBIN/mumemto" -l "$L" -o "$P/mumemto_L$L" "$R" "$Q"; }
    cut -f3- "$P/mummer_mum_L$L.set" | LC_ALL=C sort -u > "$P/mummer_mum_L$L.set3"
    if [ -s "$P/mumemto_L$L.mums" ]; then
      norm_mumemto "$P/mumemto_L$L.mums" "$P/mumemto_L$L.set"
      check "$pair" "$L" "mumemto(+,+) vs mummer -mum" "$P/mumemto_L$L.set" "$P/mummer_mum_L$L.set3" info
    else
      cached || mark_dnf "$TSV" "DNF (timeout $MUMEMTO_TIMEOUT s)"
      printf "%s\t%s\tmumemto(+,+) vs mummer -mum\tNA\t%s\tNA\t%s\tDNF\n" "$pair" "$L" \
        "$(wc -l < "$P/mummer_mum_L$L.set3")" "$(wc -l < "$P/mummer_mum_L$L.set3")" | tee -a "$CHK"
    fi

    # SMEMs from prebuilt FM-indexes (index time above).
    cached || measure "$TSV" "$pair" "$L" "bwa fastmap" "SMEM (fwd+rc)" 1 -- "$P/bwa_fastmap_L$L.txt" \
      "$ENVBIN/bwa" fastmap -l "$L" "$IDX/$pair.bwa" "$Q"
    for T in $THREADS; do
      cached || measure "$TSV" "$pair" "$L" "ropebwt3 mem" "SMEM (fwd+rc, no pos)" "$T" -- "$P/ropebwt3_L${L}_T$T.txt" \
        "$ENVBIN/ropebwt3" mem -t "$T" -l "$L" "$IDX/$pair.fmd" "$Q"
    done
  done
  # Raw outputs are large on chr1; the sorted sets are enough to re-check (SPOROS_ONLY=1 checks against them).
  rm -f "$P"/*.txt
done
echo "run_mems.sh done"
